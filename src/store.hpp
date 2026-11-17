#pragma once
#include "model.hpp"
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <queue>
#include <unordered_set>
#include <unistd.h>

namespace gc {
// Scratch ownership and bounded fan-in follow gwaslocusfinder's ExternalSort.
// See THIRD_PARTY.md for its MIT notice.
class ScratchDirectory {
    fs::path path_;
    size_t serial_ = 0;
public:
    explicit ScratchDirectory(std::string directory) {
        if (directory.empty()) {
            const char* env = std::getenv("TMPDIR");
            directory = env && *env ? env : fs::temp_directory_path().string();
        }
        std::string pattern = (fs::absolute(directory)/"genecontext.XXXXXX").string();
        require(mkdtemp(pattern.data()) != nullptr,
                "Cannot create temporary directory in " + directory + ": " + std::strerror(errno));
        path_ = pattern;
    }
    ScratchDirectory(const ScratchDirectory&) = delete;
    ScratchDirectory& operator=(const ScratchDirectory&) = delete;
    ~ScratchDirectory() { std::error_code ec; fs::remove_all(path_,ec); }
    std::string file(const std::string& suffix = ".tsv") {
        return (path_/(std::to_string(serial_++)+suffix)).string();
    }
};
class FeatureSort {
    ScratchDirectory& scratch_;
    size_t limit_, bytes_ = 0, spills_ = 0;
    std::vector<Feature> buffer_;
    std::vector<std::string> runs_;
    static void write(std::ostream& out, const Feature& feature) {
        out << feature.line();
        require(bool(out), "Cannot write temporary sort run (check free disk space)");
    }
    void spill() {
        if (buffer_.empty()) return;
        std::sort(buffer_.begin(),buffer_.end());
        std::string path = scratch_.file();
        std::ofstream out(path,std::ios::binary);
        require(bool(out), "Cannot open temporary sort run " + path);
        for (const auto& feature : buffer_) write(out,feature);
        out.close();
        require(bool(out), "Cannot close temporary sort run " + path);
        runs_.push_back(path);
        ++spills_;
        std::vector<Feature>().swap(buffer_);
        bytes_ = 0;
    }
    static void merge(const std::vector<std::string>& paths, const std::function<void(const Feature&)>& visit) {
        std::vector<std::ifstream> files(paths.size());
        std::vector<Feature> heads(paths.size());
        auto later = [&](size_t a, size_t b) { return heads[b] < heads[a]; };
        std::priority_queue<size_t,std::vector<size_t>,decltype(later)> pending(later);
        auto next = [&](size_t i) {
            std::string line;
            if (!std::getline(files[i],line)) {
                require(files[i].eof() && !files[i].bad(), "Cannot read temporary sort run " + paths[i]);
                return;
            }
            heads[i] = parse_table(line);
            pending.push(i);
        };
        for (size_t i = 0; i < paths.size(); ++i) {
            files[i].open(paths[i],std::ios::binary);
            require(bool(files[i]), "Cannot open temporary sort run " + paths[i]);
            next(i);
        }
        while (!pending.empty()) {
            size_t i = pending.top();
            pending.pop();
            visit(heads[i]);
            next(i);
        }
    }
public:
    FeatureSort(ScratchDirectory& scratch, size_t bytes) : scratch_(scratch), limit_(bytes) {}
    void add(Feature feature) {
        bytes_ += feature.memory_bytes()+sizeof(Feature); // allow for vector capacity growth
        buffer_.push_back(std::move(feature));
        if (bytes_ >= limit_) spill();
    }
    size_t spills() const { return spills_; }
    void finish(const std::function<void(const Feature&)>& visit) {
        if (runs_.empty()) {
            std::sort(buffer_.begin(),buffer_.end());
            for (const auto& feature : buffer_) visit(feature);
            std::vector<Feature>().swap(buffer_);
            bytes_ = 0;
            return;
        }
        spill();
        constexpr size_t fan_in = 32;
        while (runs_.size() > fan_in) {
            std::vector<std::string> next;
            for (size_t i = 0; i < runs_.size(); i += fan_in) {
                std::vector<std::string> group(runs_.begin()+i,runs_.begin()+std::min(runs_.size(),i+fan_in));
                std::string path = scratch_.file();
                std::ofstream out(path,std::ios::binary);
                require(bool(out), "Cannot open temporary merged run " + path);
                merge(group,[&](const Feature& feature) { write(out,feature); });
                out.close();
                require(bool(out), "Cannot close temporary merged run " + path);
                next.push_back(path);
                for (const auto& old : group) fs::remove(old);
            }
            runs_ = std::move(next);
        }
        merge(runs_,visit);
        for (const auto& path : runs_) fs::remove(path);
        runs_.clear();
    }
};
struct TbxCloser { void operator()(tbx_t* p) const { if (p) tbx_destroy(p); } };
using Tbx = std::unique_ptr<tbx_t,TbxCloser>;
struct IteratorCloser { void operator()(hts_itr_t* p) const { if (p) hts_itr_destroy(p); } };
using Iterator = std::unique_ptr<hts_itr_t,IteratorCloser>;
class BgzfWriter {
    BGZF* file_;
public:
    explicit BgzfWriter(const std::string& path) : file_(bgzf_open(path.c_str(),"w")) {
        require(file_ != nullptr, "Cannot open temporary BGZF file " + path);
    }
    BgzfWriter(const BgzfWriter&) = delete;
    BgzfWriter& operator=(const BgzfWriter&) = delete;
    ~BgzfWriter() { if (file_) bgzf_close(file_); }
    void write(const std::string& line) {
        require(bgzf_write(file_,line.data(),line.size()) == static_cast<ssize_t>(line.size()),
                "Cannot write temporary BGZF file (check free disk space)");
    }
    void close() {
        BGZF* file = file_;
        file_ = nullptr;
        require(bgzf_close(file) == 0, "Cannot close temporary BGZF file (check free disk space)");
    }
};
struct StoreOptions {
    std::string annotation, temp_dir;
    size_t sort_bytes = 64*1024*1024;
    GeneFilter filter;
    FeatureSelection selection;
    bool quiet = false;
};
class GeneStore {
    StoreOptions options_;
    std::unique_ptr<ScratchDirectory> scratch_; // destroyed after file and index
    HtsFile file_;
    Tbx index_;
    std::string path_;
    bool prepared_ = false, empty_ = false;
    std::unordered_set<std::string> source_contigs_;
    void prepare() {
        prepared_ = true;
        scratch_ = std::make_unique<ScratchDirectory>(options_.temp_dir);
        FeatureSort sorter(*scratch_,options_.sort_bytes);
        KString line;
        int status;
        size_t line_number = 0, gene_records = 0, accepted = 0;
        while ((status = hts_getline(file_.get(),'\n',&line.value)) >= 0) {
            ++line_number;
            if (fasta_directive(line.view())) break;
            try {
                std::string contig;
                auto feature = parse_feature(line.view(),options_.filter,options_.selection,&contig);
                if (!contig.empty()) source_contigs_.insert(contig);
                if (!feature) continue;
                if (feature->kind == FeatureKind::Gene) {
                    ++gene_records;
                    if (!options_.filter.accepts(feature->region)) continue;
                    ++accepted;
                }
                sorter.add(std::move(*feature));
            } catch (const std::exception& error) {
                throw std::runtime_error(options_.annotation+":"+std::to_string(line_number)+": "+error.what());
            }
        }
        require(status >= -1, "Cannot read annotation " + options_.annotation);
        require(hts_close(file_.release()) == 0, "Cannot close annotation " + options_.annotation);
        require(gene_records > 0, "No gene features found; check --features and the annotation format");
        path_ = scratch_->file(".tsv.gz");
        BgzfWriter out(path_);
        std::optional<Feature> previous;
        size_t written = 0, child_records = 0;
        sorter.finish([&](const Feature& feature) {
            if (previous && *previous == feature) return;
            out.write(feature.line());
            previous = feature;
            if (feature.kind == FeatureKind::Gene) ++written;
            else ++child_records;
        });
        out.close();
        if (!accepted) {
            empty_ = true;
            std::fprintf(stderr,"genecontext: warning: no genes pass the biotype filters\n");
        } else {
            tbx_conf_t config{TBX_GENERIC,1,2,3,'#',0};
            require(tbx_index_build(path_.c_str(),14,&config) == 0, "Cannot index temporary gene table");
            index_.reset(tbx_index_load(path_.c_str()));
            require(bool(index_), "Cannot load temporary gene index");
            file_.reset(hts_open(path_.c_str(),"r"));
            require(bool(file_), "Cannot open temporary gene table");
        }
        if (!options_.quiet)
            std::fprintf(stderr,"genecontext: prepared %zu genes from %zu gene records (%zu sort runs); %zu transcript/exon/CDS records\n",
                         written,gene_records,sorter.spills(),child_records);
    }
    bool load_existing_index() {
        bool saw_index = false;
        if (hts_get_format(file_.get())->compression == bgzf) {
            for (const auto* suffix : {".csi",".tbi"}) {
                std::string idx_path = path_+suffix;
                if (!fs::exists(idx_path)) continue;
                saw_index = true;
                Tbx candidate(tbx_index_load3(path_.c_str(),idx_path.c_str(),HTS_IDX_SILENT_FAIL));
                if (!candidate) {
                    std::fprintf(stderr,"genecontext: cannot load %s; trying another index or preparing a replacement\n",idx_path.c_str());
                    continue;
                }
                const auto& c = candidate->conf;
                if (c.preset != TBX_GENERIC || c.sc != 1 || c.bc != 4 || c.ec != 5 || c.meta_char != '#' || c.line_skip != 0) {
                    std::fprintf(stderr,"genecontext: incompatible index %s (expected 1-based columns 1,4,5); preparing if no usable index remains\n",idx_path.c_str());
                    continue;
                }
                if (fs::last_write_time(idx_path) < fs::last_write_time(path_)) {
                    std::fprintf(stderr,"genecontext: index %s is older than its annotation; preparing if no usable index remains\n",idx_path.c_str());
                    continue;
                }
                index_ = std::move(candidate);
                if (!options_.quiet) std::fprintf(stderr,"genecontext: using existing index %s\n",idx_path.c_str());
                return true;
            }
        } else if (fs::exists(path_+".csi") || fs::exists(path_+".tbi")) {
            saw_index = true;
            std::fprintf(stderr,"genecontext: annotation is not BGZF; existing index cannot be used\n");
        }
        if (!options_.quiet)
            std::fprintf(stderr,"genecontext: %s; preparing temporary gene table\n",saw_index ? "no usable index" : "no index");
        return false;
    }

    std::vector<Feature> read_features(const std::string& contig, int64_t start, int64_t end) {
        std::vector<Feature> result;
        if (empty_) return result;
        int tid = tbx_name2id(index_.get(),contig.c_str());
        if (tid < 0) return result;
        Iterator it(tbx_itr_queryi(index_.get(),tid,start-1,end));
        require(bool(it), "Cannot query annotation at " + contig+":"+std::to_string(start));
        KString line;
        int status;
        while ((status = tbx_itr_next(file_.get(),index_.get(),it.get(),&line.value)) >= 0) {
            try {
                auto feature = prepared_ ? std::optional<Feature>(parse_table(line.view())) :
                    parse_feature(line.view(),options_.filter,options_.selection);
                if (feature && feature->region.start <= end && feature->region.end >= start)
                    result.push_back(std::move(*feature));
            } catch (const std::exception& error) {
                throw std::runtime_error(path_+" near "+contig+":"+std::to_string(start)+": "+error.what());
            }
        }
        require(status == -1, "Indexed annotation read failed: " + path_);
        sort_unique(result);
        return result;
    }

public:
    explicit GeneStore(StoreOptions options) : options_(std::move(options)), path_(options_.annotation) {
        require(path_ != "-" && fs::is_regular_file(path_),
                "Annotation must be a local, rereadable GTF/GFF3 file: " + path_);
        file_.reset(hts_open(path_.c_str(),"r"));
        require(bool(file_), "Cannot open annotation " + path_);
        if (!load_existing_index()) prepare();
    }
    bool prepared() const { return prepared_; }
    bool has_contig(const std::string& contig) const {
        return prepared_ ? source_contigs_.count(contig) : tbx_name2id(index_.get(),contig.c_str()) >= 0;
    }
    GeneModel query(const std::string& contig, int64_t start, int64_t end) {
        auto features = read_features(contig,start,end);
        if (options_.selection.transcripts) {
            // Resolve transcript models across each candidate gene's full span.
            // Otherwise a narrow query can hide a conflicting transcript part
            // or omit the exons flanking an intron. Eligibility still uses the
            // original window below, including for genes found during expansion.
            int64_t full_start = start, full_end = end;
            for (const auto& feature : features) {
                const auto& gene = feature.region;
                if (feature.kind != FeatureKind::Gene || !options_.filter.accepts(gene)) continue;
                full_start = std::min(full_start,gene.start);
                full_end = std::max(full_end,gene.end);
            }
            for (const auto& [extra_start, extra_end] : {std::make_pair(full_start,start-1),std::make_pair(end+1,full_end)}) {
                if (extra_start > extra_end) continue;
                auto extra = read_features(contig,extra_start,extra_end);
                features.insert(features.end(),std::make_move_iterator(extra.begin()),std::make_move_iterator(extra.end()));
            }
            sort_unique(features);
        }
        return resolve_model(features,options_.filter,start,end);
    }
};
}
