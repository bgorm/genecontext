#pragma once
#include "common.hpp"
#include <optional>
#include <set>
#include <tuple>
#include <unordered_map>
#include <map>

namespace gc {
struct Gene {
    std::string chrom, id, name, biotype;
    int64_t start = 0, end = 0; // 1-based inclusive, as in GTF/GFF3
    char strand = '.';
    const std::string& label() const { return name.empty() ? id : name; }
    auto key() const { return std::tie(chrom,start,end,id,name,strand,biotype); }
    bool operator<(const Gene& other) const { return key() < other.key(); }
    bool operator==(const Gene& other) const { return key() == other.key(); }
    size_t memory_bytes() const {
        return sizeof(*this)+chrom.capacity()+id.capacity()+name.capacity()+biotype.capacity();
    }
    std::optional<int64_t> tss() const {
        if (strand == '+') return start;
        if (strand == '-') return end;
        return std::nullopt;
    }
    // Gene coordinate minus POS: negative to the left, positive to the right.
    int64_t signed_body_distance(int64_t pos) const {
        if (pos < start) return start - pos;
        if (pos > end) return end - pos;
        return 0;
    }
    std::string line() const {
        return chrom+'\t'+std::to_string(start)+'\t'+std::to_string(end)+'\t'+
            escape(id)+'\t'+escape(name)+'\t'+strand+'\t'+escape(biotype)+'\n';
    }
};
enum class Format { Auto, Gtf, Gff3 };
struct GeneFilter {
    Format format = Format::Auto;
    std::set<std::string> features{"gene","pseudogene","ncRNA_gene"};
    std::set<std::string> transcript_features{"transcript","mRNA","ncRNA","lnc_RNA","lncRNA",
        "tRNA","rRNA","miRNA","pre_miRNA","snRNA","snoRNA","scRNA","primary_transcript"};
    std::set<std::string> include, exclude;
    bool accepts(const Gene& gene) const {
        return (include.empty() || include.count(gene.biotype)) && !exclude.count(gene.biotype);
    }
};
using Attributes = std::unordered_map<std::string,std::vector<std::string>>;
inline void attribute(Attributes& attrs, std::string key, std::string value) {
    // Repeated GTF tags must survive: Ensembl_canonical is often not first.
    attrs[std::move(key)].push_back(std::move(value));
}
inline Attributes parse_gtf(std::string_view text) {
    Attributes result;
    size_t p = 0;
    while (p < text.size()) {
        while (p < text.size() && (text[p] == ' ' || text[p] == '\t' || text[p] == ';')) ++p;
        if (p == text.size()) break;
        size_t begin = p;
        while (p < text.size() && text[p] != ' ' && text[p] != '\t' && text[p] != ';') ++p;
        std::string key(text.substr(begin,p-begin));
        while (p < text.size() && (text[p] == ' ' || text[p] == '\t')) ++p;
        require(p < text.size() && text[p] != ';', "Missing GTF attribute value: " + key);
        std::string value;
        if (text[p] == '"') {
            ++p;
            bool closed = false;
            while (p < text.size()) {
                char c = text[p++];
                if (c == '"') { closed = true; break; }
                if (c == '\\') {
                    require(p < text.size(), "Incomplete GTF escape");
                    c = text[p++];
                }
                require(c != '\0', "NUL in GTF attribute");
                value += c;
            }
            require(closed, "Unterminated GTF attribute: " + key);
        } else {
            // GENCODE mixes quoted strings with bare scalars, e.g. level 2
            // and exon_number 1. A bare value is exactly one token.
            while (p < text.size() && text[p] != ';' && text[p] != ' ' && text[p] != '\t') {
                const char c = text[p++];
                require(c != '\0' && c != '"' && c != '\\', "Invalid unquoted GTF attribute: " + key);
                value += c;
            }
        }
        while (p < text.size() && (text[p] == ' ' || text[p] == '\t')) ++p;
        require(p == text.size() || text[p] == ';', "Missing GTF attribute separator after " + key);
        attribute(result,std::move(key),std::move(value));
    }
    return result;
}
inline Attributes parse_gff3(std::string_view text) {
    Attributes result;
    for (const auto& part : split(text,';')) {
        if (trim(part).empty()) continue;
        auto eq = part.find('=');
        require(eq != part.npos && eq != 0, "Expected key=value GFF3 attribute: " + part);
        const std::string key = trim(std::string_view(part).substr(0,eq));
        auto value = std::string_view(part).substr(eq+1);
        if (key == "Parent" || key == "tag") {
            // Split before decoding: an encoded comma is part of one ID/tag.
            for (const auto& item : split(value,',')) attribute(result,key,unescape(item));
        } else attribute(result,key,unescape(value));
    }
    return result;
}
inline std::string pick(const Attributes& attrs, std::initializer_list<const char*> keys) {
    for (const auto* key : keys) {
        auto found = attrs.find(key);
        if (found != attrs.end() && !found->second.empty() && !found->second.front().empty() && found->second.front() != ".")
            return found->second.front();
    }
    return {};
}
inline bool fasta_directive(std::string_view line) {
    return line.substr(0,7) == "##FASTA" &&
        (line.size() == 7 || line[7] == ' ' || line[7] == '\t' || line[7] == '\r');
}
enum class FeatureKind { Gene, Transcript, Exon, Cds };
struct FeatureSelection { bool transcripts = false, exons = false, cds = false; };
struct Feature {
    Gene region;
    FeatureKind kind = FeatureKind::Gene;
    std::string source_id, gene_id, transcript_id;
    std::vector<std::string> parents, tags;
    auto key() const { return std::tuple_cat(region.key(),std::tie(kind,source_id,gene_id,transcript_id,parents,tags)); }
    bool operator<(const Feature& other) const { return key() < other.key(); }
    bool operator==(const Feature& other) const { return key() == other.key(); }
    size_t memory_bytes() const {
        size_t bytes = sizeof(*this)+region.memory_bytes()-sizeof(Gene)+source_id.capacity()+gene_id.capacity()+transcript_id.capacity();
        for (const auto* values : {&parents,&tags}) {
            bytes += values->capacity()*sizeof(std::string);
            for (const auto& value : *values) bytes += value.capacity();
        }
        return bytes;
    }
    std::string line() const {
        std::string result = region.line(); result.pop_back();
        std::vector<std::string> encoded_parents, encoded_tags;
        for (const auto& parent : parents) encoded_parents.push_back(escape(parent));
        for (const auto& tag : tags) encoded_tags.push_back(escape(tag));
        return result+'\t'+std::to_string(static_cast<int>(kind))+'\t'+escape(source_id)+'\t'+
            escape(gene_id)+'\t'+escape(transcript_id)+'\t'+join(encoded_parents)+'\t'+join(encoded_tags)+'\n';
    }
};
// One parser and one internal representation serve indexed and prepared inputs.
inline std::optional<Feature> parse_feature(std::string_view line, const GeneFilter& filter,
                                          const FeatureSelection& selection, std::string* contig = nullptr) {
    if (line.empty() || line[0] == '#' || trim(line).empty()) return std::nullopt;
    if (line.back() == '\r') line.remove_suffix(1);
    auto columns = split(line,'\t');
    require(columns.size() == 9, "Expected nine tab-separated GTF/GFF3 columns");
    require(!columns[0].empty() && columns[0].find_first_of(" \r\n") == columns[0].npos, "Invalid annotation contig");
    // Seqids are kept verbatim; explicit contig mapping avoids guessing aliases.
    if (contig) *contig = columns[0];
    Feature feature;
    if (filter.features.count(columns[2])) feature.kind = FeatureKind::Gene;
    else if (selection.transcripts && filter.transcript_features.count(columns[2])) feature.kind = FeatureKind::Transcript;
    else if (selection.exons && columns[2] == "exon") feature.kind = FeatureKind::Exon;
    else if (selection.cds && columns[2] == "CDS") feature.kind = FeatureKind::Cds;
    else return std::nullopt;
    auto& region = feature.region;
    region.chrom = columns[0];
    region.start = integer(columns[3],1,HTS_POS_MAX-1,"gene start");
    region.end = integer(columns[4],region.start,HTS_POS_MAX-1,"gene end");
    require(columns[6].size() == 1 && std::string("+-.?").find(columns[6][0]) != std::string::npos, "Invalid gene strand");
    region.strand = columns[6][0] == '?' ? '.' : columns[6][0];
    Format format = filter.format;
    if (format == Format::Auto) {
        columns[8] = trim(columns[8]);
        auto p = columns[8].find_first_of("= \t");
        require(p != columns[8].npos, "Cannot detect GTF/GFF3 attributes");
        format = columns[8][p] == '=' ? Format::Gff3 : Format::Gtf;
    }
    auto attrs = format == Format::Gff3 ? parse_gff3(columns[8]) : parse_gtf(columns[8]);
    if (feature.kind == FeatureKind::Gene) region.id = pick(attrs,{"gene_id","ID","locus_tag"});
    else if (feature.kind == FeatureKind::Transcript) region.id = pick(attrs,{"transcript_id","ID"});
    else region.id = pick(attrs,{"exon_id","ID"});
    region.name = pick(attrs,{"gene_name","Name","gene"});
    region.biotype = pick(attrs,{"gene_type","gene_biotype","biotype"});
    if (feature.kind == FeatureKind::Gene) require(!region.id.empty(), "Gene feature lacks gene_id, ID or locus_tag");
    if (feature.kind == FeatureKind::Transcript) require(!region.id.empty(), "Transcript feature lacks transcript_id or ID");
    feature.source_id = pick(attrs,{"ID"});
    feature.gene_id = pick(attrs,{"gene_id"});
    feature.transcript_id = pick(attrs,{"transcript_id"});
    if (attrs.count("Parent")) feature.parents = attrs.at("Parent");
    if (attrs.count("tag")) feature.tags = attrs.at("tag");
    sort_unique(feature.parents);
    sort_unique(feature.tags);
    return feature;
}
inline Feature parse_table(std::string_view line) {
    auto columns = split(line,'\t');
    require(columns.size() == 13 && columns[5].size() == 1, "Invalid temporary annotation table");
    Feature feature;
    auto& region = feature.region;
    region.chrom = columns[0];
    region.start = integer(columns[1],1,HTS_POS_MAX-1,"temporary gene start");
    region.end = integer(columns[2],region.start,HTS_POS_MAX-1,"temporary gene end");
    region.id = unescape(columns[3]);
    region.name = unescape(columns[4]);
    region.strand = columns[5][0];
    region.biotype = unescape(columns[6]);
    feature.kind = static_cast<FeatureKind>(integer(columns[7],0,3,"temporary feature type"));
    feature.source_id = unescape(columns[8]);
    feature.gene_id = unescape(columns[9]);
    feature.transcript_id = unescape(columns[10]);
    if (!columns[11].empty()) for (const auto& value : split(columns[11],',')) feature.parents.push_back(unescape(value));
    if (!columns[12].empty()) for (const auto& value : split(columns[12],',')) feature.tags.push_back(unescape(value));
    return feature;
}
}
