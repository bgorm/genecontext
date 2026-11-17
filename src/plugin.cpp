#include "options.hpp"
#include <htslib/vcf.h>
#include <sstream>

#ifndef GC_BCFTOOLS_VERSION
#error Build with the target bcftools version
#endif
#ifndef GC_HTSLIB_VERSION
#error Build with the target HTSlib version
#endif
#define GC_VERSION "0.2.3"
#define GC_EXPORT extern "C" __attribute__((visibility("default")))
extern "C" const char* bcftools_version(void);

namespace gc {
static const char help[] =
    "About: Annotate variant gene neighborhoods from GTF/GFF3 in a bcftools stream.\n"
    "(genecontext version " GC_VERSION ")\n\n"
    "Usage: bcftools +genecontext [General Options] <in.vcf|in.bcf|-> -- [Plugin Options]\n\n"
    "   bcftools handles -i/-e, -r/-R, -t/-T, -o, -O, --threads and -W before --.\n"
    "   All new annotations are INFO fields and use the variant's 1-based POS.\n\n"
    "Plugin options (after --):\n"
    "   -a, --annotation FILE      Gene annotation: plain/gzip/BGZF GTF or GFF3 (required)\n"
    "       --emit LIST            Comma-separated output groups [nearest,distance,nearest-tss,tss-distance]\n"
    "                              nearest, distance, nearest-tss, tss-distance, neighbors,\n"
    "                              neighbors-tss, map, overlap, ids, biotype, strand,\n"
    "                              flanks, count, overlap-count, tss-details, promoter,\n"
    "                              exon, cds, intron, feature-details, relations; or all\n"
    "       --tss-source gene|canonical  Gene boundaries or tagged transcript TSSs [gene]\n"
    "       --canonical-tag STR    Transcript tag used for canonical selection [Ensembl_canonical]\n"
    "       --promoter-upstream INT   Promoter bp upstream of the selected TSS [2000]\n"
    "       --promoter-downstream INT Promoter bp downstream of the selected TSS [200]\n"
    "       --feature-transcripts all|canonical  Transcripts used for exon/CDS/intron output [all]\n"
    "       --max-distance INT     Maximum nearest-gene/TSS/flank distance in bp [1000000]\n"
    "       --window INT           Neighborhood radius on each side of POS in bp [250000]\n"
    "       --dot-bp INT           Base pairs per dot in the map, rounded up [10000]\n"
    "       --max-dots INT         Maximum dots per map gap; ~ marks longer gaps [80]\n"
    "       --label name|id        Gene labels; missing names always fall back to IDs [name]\n"
    "       --protein-coding       Shortcut for --biotypes protein_coding\n"
    "       --biotypes LIST        Include exact comma-separated gene biotypes [all]\n"
    "       --exclude-biotypes LIST  Exclude exact comma-separated gene biotypes\n"
    "       --features LIST        Gene-level feature types [gene,pseudogene,ncRNA_gene]\n"
    "       --transcript-features LIST  Transcript feature types [transcript,mRNA and common RNA types]\n"
    "       --annotation-format auto|gtf|gff3  Attribute syntax [auto]\n"
    "       --chr-map FILE         Two columns: VCF contig, annotation contig [exact names]\n"
    "       --missing-contig MODE  warn, error or ignore [warn]\n"
    "       --prefix STR           Prefix for emitted INFO tags [GC_]\n"
    "       --overwrite            Replace selected existing tags of compatible Number/Type\n"
    "       --sort-mem INT         Approximate sort buffer in MiB for preparation [64]\n"
    "       --temp-dir DIR         Scratch directory parent [TMPDIR or system temporary directory]\n"
    "       --quiet                Suppress preparation and completion messages\n"
    "   -h, --help                 Show this help\n\n"
    "Definitions:\n"
    "   Gene body includes both endpoints. Distance is zero inside the body. All ties survive.\n"
    "   TSS uses the gene or tagged transcript start on + and end on -. Unknown strands have no TSS.\n"
    "   Canonical mode has no fallback. All transcripts bearing the chosen tag are eligible.\n"
    "   Promoter offsets are strand-aware; introns are gaps between merged exons of one transcript.\n"
    "   Transcript features are required for canonical TSSs; exon/CDS records supply feature overlaps.\n"
    "   Neighborhood distances are negative to the genomic left, positive to the right.\n"
    "   Map anchors are the nearest body boundary or POS for [overlapping|genes]; [] is intergenic.\n"
    "   Dots measure gaps between anchors and window edges. Tied anchors use |.\n"
    "   Existing compatible CSI/TBI indexes are reused; otherwise a temporary table is prepared.\n"
    "   Missing nearest matches remain missing; searches do not extend beyond max-distance.\n\n"
    "Examples:\n"
    "   bcftools +genecontext index.bcf -Ob -o genes.bcf -W -- -a genes.gtf.gz\n"
    "   bcftools view -Ou input.bcf | bcftools +genecontext -Oz -o genes.vcf.gz -- \\\n"
    "       -a genes.gff3.gz --emit all --protein-coding --window 100000 --dot-bp 10000\n"
    "   bcftools +genecontext index.bcf -- -a full.gtf.gz --tss-source canonical \\\n"
    "       --emit nearest-tss,tss-distance,tss-details,promoter,exon,cds,intron,feature-details\n";

std::string quoted_header_value(const std::string& text) {
    std::string result = "\"";
    for (const char c : text) {
        if (c == '\\' || c == '"') result += '\\';
        result += (c == '\n' || c == '\r' || c == '\t') ? ' ' : c;
    }
    return result + '"';
}

std::string provenance_header(const Options& options, bool prepared) {
    const auto& annotation = options.annotation;
    const auto& filter = options.store.filter;
    std::ostringstream header;
    header.imbue(std::locale::classic());
    header << "##genecontext=<Version=" << quoted_header_value(GC_VERSION)
           << ",Annotation=" << quoted_header_value(fs::absolute(options.store.annotation).string())
           << ",Source=" << quoted_header_value(prepared ? "prepared" : "indexed")
           << ",Coordinate=\"POS\""
           << ",TSS=" << quoted_header_value(annotation.canonical_tss ? "canonical_transcript" : "gene_boundary")
           << ",CanonicalTag=" << quoted_header_value(annotation.canonical_tag)
           << ",FeatureTranscripts=" << quoted_header_value(annotation.canonical_features ? "canonical" : "all")
           << ",PromoterUpstream=" << annotation.promoter_upstream
           << ",PromoterDownstream=" << annotation.promoter_downstream
           << ",MaxDistance=" << annotation.max_distance
           << ",Window=" << annotation.window
           << ",DotBP=" << annotation.dot_bp
           << ",MaxDots=" << annotation.max_dots
           << ",Biotypes=" << quoted_header_value(filter.include.empty() ? "all" : join(filter.include))
           << ",ExcludeBiotypes=" << quoted_header_value(join(filter.exclude))
           << ",Prefix=" << quoted_header_value(annotation.prefix)
           << ",Labels=" << quoted_header_value(annotation.label_ids ? "id" : "name")
           << ",Emit=" << quoted_header_value(join(annotation.emit))
           << ",GeneFeatures=" << quoted_header_value(join(filter.features))
           << ",TranscriptFeatures=" << quoted_header_value(join(filter.transcript_features))
           << ",ContigMap=" << quoted_header_value(options.chr_map) << '>';
    return header.str();
}

std::unordered_map<std::string,std::string> read_contig_map(const std::string& path) {
    std::unordered_map<std::string,std::string> mapping;
    if (path.empty()) return mapping;
    std::ifstream input(path);
    require(bool(input), "Cannot open contig map " + path);
    for (std::string line; std::getline(input, line);) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        std::istringstream columns(line);
        std::string from, to, extra;
        require(bool(columns >> from >> to) && !(columns >> extra),
                "Contig map must have two columns: VCF contig, annotation contig");
        require(mapping.emplace(from, to).second, "Duplicate VCF contig in --chr-map: " + from);
    }
    require(input.eof() && !input.bad(), "Cannot read contig map " + path);
    return mapping;
}

class Context {
    struct OutputField {
        const Field* definition;
        std::string tag;
    };

    Options options_;
    bcf_hdr_t* header_; // owned by bcftools
    std::unique_ptr<GeneStore> store_;
    std::vector<OutputField> output_fields_;
    std::unordered_map<std::string,std::string> chr_map_;
    std::set<std::string> missing_contigs_;
    uint64_t records_ = 0;
    int last_rid_ = -1;
    int64_t last_pos_ = -1;
    Annotations cached_;

    void add_info_headers() {
        for (const auto& field : fields) {
            if (!options_.annotation.emit.count(field.option)) continue;
            const std::string tag = options_.annotation.prefix + field.suffix;
            const int id = bcf_hdr_id2int(header_, BCF_DT_ID, tag.c_str());
            if (id >= 0 && bcf_hdr_idinfo_exists(header_, BCF_HL_INFO, id)) {
                require(options_.overwrite, "INFO/" + tag + " already exists; use --overwrite or another --prefix");
                const bool compatible_type = bcf_hdr_id2type(header_, BCF_HL_INFO, id) ==
                                             (field.is_string() ? BCF_HT_STR : BCF_HT_INT);
                const bool compatible_number = bcf_hdr_id2length(header_, BCF_HL_INFO, id) ==
                                               (field.is_variable() ? BCF_VL_VAR : BCF_VL_FIXED) &&
                                               (field.is_variable() || bcf_hdr_id2number(header_, BCF_HL_INFO, id) == 1);
                require(compatible_type && compatible_number, "Incompatible Number/Type for INFO/" + tag);
                bcf_hdr_remove(header_, BCF_HL_INFO, tag.c_str());
            }
            const std::string line = "##INFO=<ID=" + tag + ",Number=" + field.number + ",Type=" + field.type +
                                     ",Description=\"" + field.description + "\">";
            require(bcf_hdr_append(header_, line.c_str()) == 0, "Cannot add INFO header " + tag);
            output_fields_.push_back({&field, tag});
        }
        require(bcf_hdr_sync(header_) == 0, "Cannot synchronize INFO headers");
    }

    Annotations query_annotations(int rid, int64_t zero_based_pos) {
        const char* raw = bcf_hdr_id2name(header_, rid);
        require(raw != nullptr, "Variant contig missing from header");
        std::string chrom = raw;
        const auto mapped = chr_map_.find(chrom);
        if (mapped != chr_map_.end()) chrom = mapped->second;

        if (!store_->has_contig(chrom)) {
            require(options_.missing_contig != "error",
                    "Annotation has no contig " + chrom + "; use --chr-map for explicit aliases");
            if (missing_contigs_.insert(chrom).second && options_.missing_contig == "warn")
                std::fprintf(stderr, "genecontext: warning: annotation has no contig %s; leaving INFO missing (check --chr-map)\n",
                             chrom.c_str());
            return {};
        }

        const int64_t pos = zero_based_pos + 1;
        const int64_t radius = options_.annotation.radius();
        const auto model = store_->query(chrom, std::max<int64_t>(1, pos - radius), pos + radius);
        return annotate(model, pos, options_.annotation);
    }

    void update_info(bcf1_t* record) const {
        for (const auto& output : output_fields_) {
            const auto& field = *output.definition;
            int status;
            if (field.is_string()) {
                const auto found = cached_.strings.find(field.suffix);
                const char* value = found == cached_.strings.end() || found->second.empty() ? nullptr : found->second.c_str();
                status = bcf_update_info_string(header_, record, output.tag.c_str(), value);
            } else {
                const auto found = cached_.integers.find(field.suffix);
                const int32_t* value = found == cached_.integers.end() ? nullptr : &found->second;
                status = bcf_update_info_int32(header_, record, output.tag.c_str(), value, value ? 1 : 0);
            }
            require(status == 0, "Cannot update INFO/" + output.tag);
        }
    }

public:
    Context(Options options, bcf_hdr_t* header) : options_(std::move(options)), header_(header) {
        add_info_headers(); // reject collisions before expensive preparation
        chr_map_ = read_contig_map(options_.chr_map);
        store_ = std::make_unique<GeneStore>(options_.store);
        const std::string provenance = provenance_header(options_, store_->prepared());
        require(bcf_hdr_append(header_, provenance.c_str()) == 0 && bcf_hdr_sync(header_) == 0,
                "Cannot add provenance header");
    }

    bcf1_t* process(bcf1_t* record) {
        require(record->rid >= 0 && record->pos >= 0 && record->pos < HTS_POS_MAX - INT32_MAX - 1,
                "Invalid variant contig or position");
        if (record->rid != last_rid_ || record->pos != last_pos_) {
            cached_ = query_annotations(record->rid, record->pos);
            last_rid_ = record->rid;
            last_pos_ = record->pos;
        }
        update_info(record);
        ++records_;
        return record;
    }

    void finish() const {
        if (!options_.store.quiet)
            std::fprintf(stderr, "genecontext: annotated %llu records\n", static_cast<unsigned long long>(records_));
    }
};
static std::unique_ptr<Context> context;
}

GC_EXPORT const char* about(void) { return "Annotate nearest genes, TSSs and gene neighborhoods from GTF/GFF3 in INFO.\n"; }
GC_EXPORT const char* usage(void) { return gc::help; }
GC_EXPORT void version(const char** bcftools, const char** htslib) {
    *bcftools = GC_BCFTOOLS_VERSION;
    *htslib = GC_HTSLIB_VERSION;
}
GC_EXPORT int init(int argc, char** argv, bcf_hdr_t*, bcf_hdr_t* out) {
    try {
        gc::require(!std::strcmp(hts_version(),GC_HTSLIB_VERSION) && !std::strcmp(bcftools_version(),GC_BCFTOOLS_VERSION),
                    "Plugin built for bcftools " GC_BCFTOOLS_VERSION " / HTSlib " GC_HTSLIB_VERSION "; rebuild for this host");
        gc::context = std::make_unique<gc::Context>(gc::parse_options(argc,argv,gc::help),out);
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr,"genecontext: %s\n",error.what());
        gc::context.reset();
        return -1;
    }
}
GC_EXPORT bcf1_t* process(bcf1_t* record) {
    try {
        gc::require(bool(gc::context),"Plugin was not initialized");
        return gc::context->process(record);
    } catch (const std::exception& error) {
        std::fprintf(stderr,"genecontext: %s\n",error.what());
        gc::context.reset();
        // The streaming API has no error return: NULL would silently discard a
        // record. A processing error must terminate bcftools with nonzero status.
        std::exit(EXIT_FAILURE);
    }
}
GC_EXPORT void destroy(void) {
    if (gc::context) gc::context->finish();
    gc::context.reset();
}
