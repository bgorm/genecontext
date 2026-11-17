#pragma once
#include "model.hpp"
#include <array>
#include <map>

namespace gc {
struct Field {
    const char *option, *suffix, *number, *type, *description;
    bool is_string() const { return std::string_view(type) == "String"; }
    bool is_variable() const { return std::string_view(number) == "."; }
};
inline const std::array<Field,28> fields{{
    {"nearest", "NEAREST", ".", "String", "Nearest gene bodies to POS within max-distance; all ties in genomic order"},
    {"distance", "DISTANCE", "1", "Integer", "Absolute distance in bp from POS to the nearest gene body; zero inside a gene"},
    {"nearest-tss", "TSS", ".", "String", "Genes with the nearest selected TSS to POS within max-distance; all ties; selection recorded in genecontext header"},
    {"tss-distance", "TSS_DISTANCE", "1", "Integer", "Absolute distance in bp from POS to the nearest selected TSS"},
    {"neighbors", "NEIGHBORS", ".", "String", "Genes within window as GENE(signed-bp-to-body); negative left and positive right in genomic coordinates"},
    {"neighbors-tss", "TSS_NEIGHBORS", ".", "String", "Genes with TSS within window as GENE(signed-bp-to-TSS); negative left and positive right in genomic coordinates"},
    {"map", "MAP", "1", "String", "Gene neighborhood at nearest body boundaries; brackets locate POS; dots are ceil(bp/dot-bp) with ~ marking capped gaps"},
    {"overlap", "OVERLAP", ".", "String", "Genes whose bodies contain POS; distinct gene IDs in genomic order"},
    {"ids", "NEAREST_ID", ".", "String", "Gene identifiers for nearest bodies; aligned with NEAREST"},
    {"ids", "TSS_ID", ".", "String", "Gene identifiers for nearest selected TSSs; aligned with TSS"},
    {"ids", "OVERLAP_ID", ".", "String", "Gene identifiers for bodies containing POS; aligned with OVERLAP"},
    {"biotype", "BIOTYPE", ".", "String", "Biotypes of nearest gene bodies; aligned with NEAREST; dot if unavailable"},
    {"strand", "STRAND", ".", "String", "Strands of nearest gene bodies; aligned with NEAREST; dot if unknown"},
    {"flanks", "LEFT", ".", "String", "Nearest nonoverlapping genes on the genomic left within max-distance; all ties"},
    {"flanks", "LEFT_DISTANCE", "1", "Integer", "Absolute bp distance from POS to the nearest gene body wholly on its genomic left"},
    {"flanks", "RIGHT", ".", "String", "Nearest nonoverlapping genes on the genomic right within max-distance; all ties"},
    {"flanks", "RIGHT_DISTANCE", "1", "Integer", "Absolute bp distance from POS to the nearest gene body wholly on its genomic right"},
    {"count", "N_GENES", "1", "Integer", "Number of distinct eligible gene IDs with body distance at most window from POS"},
    {"overlap-count", "N_OVERLAP", "1", "Integer", "Number of distinct eligible gene IDs whose bodies contain POS"},
    {"tss-details", "TSS_DETAIL", ".", "String", "Nearest TSS sites as GENE|GENE_ID|TRANSCRIPT_ID|TSS|STRAND|SOURCE; all tied sites; transcript dot for gene boundaries"},
    {"promoter", "PROMOTER", ".", "String", "Genes whose selected strand-aware TSS promoter windows contain POS"},
    {"promoter", "PROMOTER_ID", ".", "String", "Gene identifiers aligned with PROMOTER"},
    {"promoter", "PROMOTER_DETAIL", ".", "String", "Promoter hits as GENE|GENE_ID|TRANSCRIPT_ID|TSS|STRAND|OFFSET; offset is negative upstream and positive downstream in transcription direction"},
    {"exon", "EXON", ".", "String", "Genes with an annotated exon containing POS; transcript selection recorded in genecontext header"},
    {"cds", "CDS", ".", "String", "Genes with an annotated CDS containing POS; this is positional overlap not a predicted consequence"},
    {"intron", "INTRON", ".", "String", "Genes with POS in a gap between consecutive merged exons of a single transcript; transcript tails excluded"},
    {"feature-details", "FEATURE", ".", "String", "Overlapping features as GENE|GENE_ID|TRANSCRIPT_ID|FEATURE|START|END; one entry per transcript/interval; introns derived from exons"},
    {"relations", "RELATION", ".", "String", "Body-neighborhood genes as GENE(upstream/downstream/overlap/unknown_strand); relationship describes POS relative to gene strand"}
}};
struct AnnotationOptions {
    std::set<std::string> emit{"nearest","distance","nearest-tss","tss-distance"};
    int64_t max_distance = 1000000, window = 250000, dot_bp = 10000;
    size_t max_dots = 80;
    bool label_ids = false;
    bool canonical_tss = false, canonical_features = false;
    std::string canonical_tag = "Ensembl_canonical";
    int64_t promoter_upstream = 2000, promoter_downstream = 200;
    std::string prefix = "GC_";
    std::string label(const Gene& gene, bool map = false) const {
        return escape(label_ids ? gene.id : gene.label(),map);
    }
    int64_t radius() const {
        int64_t radius = 0;
        for (const auto& option : emit) {
            if (option == "neighbors" || option == "neighbors-tss" || option == "map" || option == "count" || option == "relations")
                radius = std::max(radius,window);
            else if (option == "promoter") radius = std::max({radius,promoter_upstream,promoter_downstream});
            else if (option != "overlap" && option != "overlap-count" && option != "exon" &&
                     option != "cds" && option != "intron" && option != "feature-details") radius = std::max(radius,max_distance);
        }
        return radius;
    }
    bool intron_context() const { return emit.count("intron") || emit.count("feature-details"); }
    FeatureSelection selection() const {
        FeatureSelection out;
        out.exons = emit.count("exon") || intron_context();
        out.cds = emit.count("cds") || emit.count("feature-details");
        out.transcripts = out.exons || out.cds || (canonical_tss &&
            (emit.count("nearest-tss") || emit.count("tss-distance") || emit.count("neighbors-tss") ||
             emit.count("ids") || emit.count("tss-details") || emit.count("promoter")));
        return out;
    }
};
struct Hit { const Gene* gene; int64_t distance; };
enum class GenomicSide { Either, Left, Right };
inline int64_t magnitude(int64_t value) { return value < 0 ? -value : value; }
// Multiple rows for one ID are treated as parts of one gene; the closest part
// supplies its annotation. Never join different IDs just because names coincide.
inline std::vector<Hit> body_hits(const std::vector<Gene>& genes, int64_t pos, int64_t limit,
                                  GenomicSide side = GenomicSide::Either) {
    std::unordered_map<std::string,size_t> by_id;
    std::vector<Hit> result;
    for (const auto& gene : genes) {
        const int64_t distance = gene.signed_body_distance(pos);
        if (magnitude(distance) > limit) continue;
        if (side == GenomicSide::Left && distance >= 0) continue;
        if (side == GenomicSide::Right && distance <= 0) continue;
        auto found = by_id.find(gene.id);
        if (found == by_id.end()) {
            by_id.emplace(gene.id,result.size()); result.push_back({&gene,distance});
        } else {
            auto& previous = result[found->second];
            if (magnitude(distance) < magnitude(previous.distance) ||
                (magnitude(distance) == magnitude(previous.distance) && gene < *previous.gene)) previous = {&gene,distance};
        }
    }
    std::sort(result.begin(),result.end(),[](const Hit& a, const Hit& b) {
        if (a.distance != b.distance) return a.distance < b.distance;
        return *a.gene < *b.gene;
    });
    return result;
}
inline std::vector<Hit> closest(std::vector<Hit> values) {
    if (values.empty()) return values;
    int64_t best = magnitude(values.front().distance);
    for (const auto& hit : values) best = std::min(best,magnitude(hit.distance));
    values.erase(std::remove_if(values.begin(),values.end(),[&](const Hit& hit) {
        return magnitude(hit.distance) != best;
    }),values.end());
    std::sort(values.begin(),values.end(),[](const Hit& a, const Hit& b) { return *a.gene < *b.gene; });
    return values;
}
inline std::string labels(const std::vector<Hit>& values, const AnnotationOptions& options, bool distances = false) {
    std::vector<std::string> text;
    for (const auto& hit : values) {
        std::string value = options.label(*hit.gene);
        if (distances) value += "("+std::string(hit.distance > 0 ? "+" : "")+std::to_string(hit.distance)+")";
        text.push_back(std::move(value));
    }
    return join(text);
}
inline std::string ids(const std::vector<Hit>& values) {
    std::vector<std::string> text;
    for (const auto& hit : values) text.push_back(escape(hit.gene->id));
    return join(text);
}
inline std::string dots(int64_t gap, const AnnotationOptions& options) {
    if (gap <= 0) return {};
    auto count = static_cast<uint64_t>((gap-1)/options.dot_bp+1);
    return std::string(std::min<uint64_t>(count,options.max_dots),'.')+(count > options.max_dots ? "~" : "");
}
inline std::string neighborhood_map(const std::vector<Hit>& neighbors, int64_t pos, const AnnotationOptions& options) {
    std::map<int64_t,std::vector<std::string>> anchors;
    std::vector<std::string> overlapping;
    for (const auto& hit : neighbors) {
        if (hit.distance == 0) overlapping.push_back(options.label(*hit.gene,true));
        else anchors[pos+hit.distance].push_back(options.label(*hit.gene,true));
    }
    anchors[pos].push_back("["+join(overlapping,"|")+"]");
    int64_t cursor = std::max<int64_t>(1,pos-options.window);
    std::string result;
    for (const auto& entry : anchors) {
        result += dots(entry.first-cursor,options)+join(entry.second,"|");
        cursor = entry.first;
    }
    result += dots(pos+options.window-cursor,options);
    return result;
}
struct Annotations {
    std::map<std::string,std::string> strings;
    std::map<std::string,int32_t> integers;
};
struct TssSite {
    const Gene* gene;
    std::string transcript;
    int64_t position;
    char strand;
};
inline std::vector<TssSite> tss_sites(const GeneModel& model, const AnnotationOptions& options) {
    std::vector<TssSite> result;
    if (!options.canonical_tss) {
        for (const auto& gene : model.genes) if (gene.tss()) result.push_back({&gene,"",*gene.tss(),gene.strand});
    } else {
        std::unordered_map<std::string,const Gene*> by_id;
        for (const auto& gene : model.genes) by_id.emplace(gene.id,&gene);
        for (const auto& tx : model.transcripts) {
            if (!tx.tagged(options.canonical_tag) || !tx.tss()) continue;
            auto parent = by_id.find(tx.gene_id);
            if (parent != by_id.end()) result.push_back({parent->second,tx.id,*tx.tss(),tx.strand});
        }
    }
    std::sort(result.begin(),result.end(),[](const TssSite& a, const TssSite& b) {
        return std::tuple_cat(a.gene->key(),std::tie(a.position,a.transcript,a.strand)) <
               std::tuple_cat(b.gene->key(),std::tie(b.position,b.transcript,b.strand));
    });
    result.erase(std::unique(result.begin(),result.end(),[](const TssSite& a, const TssSite& b) {
        return a.gene->id == b.gene->id && a.position == b.position && a.transcript == b.transcript && a.strand == b.strand;
    }),result.end());
    return result;
}
inline std::vector<TssSite> nearest_sites(const std::vector<TssSite>& sites, int64_t pos, int64_t limit) {
    std::vector<TssSite> result;
    int64_t best = limit;
    for (const auto& site : sites) {
        int64_t distance = magnitude(site.position-pos);
        if (distance > best) continue;
        if (distance < best) { result.clear(); best = distance; }
        result.push_back(site);
    }
    return result;
}
enum class HitOrder { GeneCoordinates, SignedDistance };
inline std::vector<Hit> site_genes(const std::vector<TssSite>& sites, int64_t pos, int64_t limit, HitOrder order = HitOrder::GeneCoordinates) {
    std::unordered_map<std::string,size_t> by_id;
    std::vector<Hit> result;
    for (const auto& site : sites) {
        int64_t distance = site.position-pos;
        if (magnitude(distance) > limit) continue;
        auto found = by_id.emplace(site.gene->id,result.size());
        if (found.second) result.push_back({site.gene,distance});
        else {
            auto& previous = result[found.first->second];
            // Sites arrive in deterministic gene/position/transcript order.
            // Keep that order when one gene has equidistant sites.
            if (magnitude(distance) < magnitude(previous.distance)) previous = {site.gene,distance};
        }
    }
    std::sort(result.begin(),result.end(),[&](const Hit& a, const Hit& b) {
        if (order == HitOrder::SignedDistance && a.distance != b.distance) return a.distance < b.distance;
        return *a.gene < *b.gene;
    });
    return result;
}
inline std::string site_detail(const TssSite& site, const AnnotationOptions& options) {
    return options.label(*site.gene)+'|'+escape(site.gene->id)+'|'+
        (site.transcript.empty() ? "." : escape(site.transcript))+'|'+std::to_string(site.position)+'|'+site.strand;
}
struct FeatureOverlap {
    const Gene* gene;
    std::string transcript, type;
    int64_t start, end;
};

// Merge overlapping/adjacent exons before looking for an internal gap.
// A position can lie in at most one intron of a given transcript.
inline std::optional<std::pair<int64_t,int64_t>> intron_at_position(
        std::vector<std::pair<int64_t,int64_t>>& exons, int64_t pos) {
    if (exons.empty()) return std::nullopt;
    std::sort(exons.begin(), exons.end());
    int64_t merged_end = exons.front().second;
    for (size_t i = 1; i < exons.size(); ++i) {
        const auto [start, end] = exons[i];
        if (merged_end + 1 < start && pos > merged_end && pos < start)
            return std::make_pair(merged_end + 1, start - 1);
        merged_end = std::max(merged_end, end);
    }
    return std::nullopt;
}

inline std::vector<FeatureOverlap> feature_overlaps(const GeneModel& model, int64_t pos,
                                                   const AnnotationOptions& options) {
    std::unordered_map<std::string,const Gene*> genes;
    for (const auto& gene : model.genes) {
        if (gene.start <= pos && gene.end >= pos) genes.emplace(gene.id, &gene);
    }

    using TranscriptKey = std::pair<std::string,std::string>; // gene ID, transcript ID
    std::set<TranscriptKey> canonical;
    if (options.canonical_features) {
        for (const auto& transcript : model.transcripts) {
            if (transcript.tagged(options.canonical_tag))
                canonical.emplace(transcript.gene_id, transcript.id);
        }
    }

    std::map<TranscriptKey,std::vector<std::pair<int64_t,int64_t>>> exons;
    std::vector<FeatureOverlap> overlaps;
    for (const auto& segment : model.segments) {
        const auto gene = genes.find(segment.gene_id);
        const TranscriptKey key{segment.gene_id, segment.transcript_id};
        if (gene == genes.end()) continue;
        if (options.canonical_features && !canonical.count(key)) continue;

        if (segment.kind == FeatureKind::Exon && !segment.transcript_id.empty() && options.intron_context())
            exons[key].emplace_back(segment.start, segment.end);

        if (pos < segment.start || pos > segment.end) continue;
        const char* type = segment.kind == FeatureKind::Exon ? "exon" : "CDS";
        overlaps.push_back({gene->second, segment.transcript_id, type, segment.start, segment.end});
    }

    for (auto& [key, intervals] : exons) {
        const auto intron = intron_at_position(intervals, pos);
        if (intron)
            overlaps.push_back({genes.at(key.first), key.second, "intron", intron->first, intron->second});
    }
    std::sort(overlaps.begin(), overlaps.end(), [](const FeatureOverlap& a, const FeatureOverlap& b) {
        return std::tuple_cat(a.gene->key(), std::tie(a.transcript, a.type, a.start, a.end)) <
               std::tuple_cat(b.gene->key(), std::tie(b.transcript, b.type, b.start, b.end));
    });
    return overlaps;
}

inline std::string feature_detail(const FeatureOverlap& overlap, const AnnotationOptions& options) {
    return options.label(*overlap.gene) + '|' + escape(overlap.gene->id) + '|' +
        (overlap.transcript.empty() ? "." : escape(overlap.transcript)) + '|' + overlap.type + '|' +
        std::to_string(overlap.start) + '|' + std::to_string(overlap.end);
}

inline void annotate_features(Annotations& out, const GeneModel& model, int64_t pos,
                              const AnnotationOptions& options) {
    const auto selection = options.selection();
    if (!selection.exons && !selection.cds) return;

    std::map<std::string,std::set<std::string>> seen;
    std::map<std::string,std::vector<std::string>> gene_labels;
    std::vector<std::string> details;
    for (const auto& overlap : feature_overlaps(model, pos, options)) {
        const std::string key = overlap.type == "exon" ? "EXON" : overlap.type == "intron" ? "INTRON" : "CDS";
        if (seen[key].insert(overlap.gene->id).second)
            gene_labels[key].push_back(options.label(*overlap.gene));
        details.push_back(feature_detail(overlap, options));
    }
    for (const auto& [key, values] : gene_labels) out.strings[key] = join(values);
    out.strings["FEATURE"] = join(details);
}

inline void annotate_gene_bodies(Annotations& out, const GeneModel& model, int64_t pos,
                                 const AnnotationOptions& options) {
    const auto nearest = closest(body_hits(model.genes, pos, options.max_distance));
    const auto overlaps = closest(body_hits(model.genes, pos, 0));
    out.strings["NEAREST"] = labels(nearest, options);
    out.strings["NEAREST_ID"] = ids(nearest);
    out.strings["OVERLAP"] = labels(overlaps, options);
    out.strings["OVERLAP_ID"] = ids(overlaps);
    if (!nearest.empty())
        out.integers["DISTANCE"] = static_cast<int32_t>(magnitude(nearest.front().distance));

    require(overlaps.size() <= INT32_MAX, "Too many genes for VCF Integer count");
    out.integers["N_OVERLAP"] = static_cast<int32_t>(overlaps.size());

    if (options.emit.count("biotype") || options.emit.count("strand")) {
        std::vector<std::string> biotypes, strands;
        for (const auto& hit : nearest) {
            biotypes.push_back(hit.gene->biotype.empty() ? "." : escape(hit.gene->biotype));
            strands.emplace_back(1, hit.gene->strand);
        }
        out.strings["BIOTYPE"] = join(biotypes);
        out.strings["STRAND"] = join(strands);
    }
    if (options.emit.count("flanks")) {
        for (const auto side : {GenomicSide::Left, GenomicSide::Right}) {
            const auto flank = closest(body_hits(model.genes, pos, options.max_distance, side));
            const std::string key = side == GenomicSide::Left ? "LEFT" : "RIGHT";
            out.strings[key] = labels(flank, options);
            if (!flank.empty())
                out.integers[key + "_DISTANCE"] = static_cast<int32_t>(magnitude(flank.front().distance));
        }
    }
}

inline const char* gene_relationship(const Hit& hit) {
    if (hit.distance == 0) return "overlap";
    if (hit.gene->strand == '.') return "unknown_strand";
    const bool variant_is_upstream = (hit.distance > 0) == (hit.gene->strand == '+');
    return variant_is_upstream ? "upstream" : "downstream";
}

inline void annotate_neighborhood(Annotations& out, const GeneModel& model, int64_t pos,
                                  const AnnotationOptions& options) {
    const auto neighbors = body_hits(model.genes, pos, options.window);
    require(neighbors.size() <= INT32_MAX, "Too many genes for VCF Integer count");
    out.integers["N_GENES"] = static_cast<int32_t>(neighbors.size());

    if (options.emit.count("neighbors")) out.strings["NEIGHBORS"] = labels(neighbors, options, true);
    if (options.emit.count("map")) out.strings["MAP"] = neighborhood_map(neighbors, pos, options);
    if (options.emit.count("relations")) {
        std::vector<std::string> relationships;
        for (const auto& hit : neighbors)
            relationships.push_back(options.label(*hit.gene) + '(' + gene_relationship(hit) + ')');
        out.strings["RELATION"] = join(relationships);
    }
}

inline void annotate_promoters(Annotations& out, const std::vector<TssSite>& sites, int64_t pos,
                               const AnnotationOptions& options) {
    std::vector<TssSite> promoters;
    std::vector<std::string> details;
    for (const auto& site : sites) {
        const int64_t offset = (pos - site.position) * (site.strand == '+' ? 1 : -1);
        if (offset < -options.promoter_upstream || offset > options.promoter_downstream) continue;
        promoters.push_back(site);
        details.push_back(site_detail(site, options) + '|' + std::to_string(offset));
    }
    const auto genes = site_genes(promoters, pos, std::max(options.promoter_upstream, options.promoter_downstream));
    out.strings["PROMOTER"] = labels(genes, options);
    out.strings["PROMOTER_ID"] = ids(genes);
    out.strings["PROMOTER_DETAIL"] = join(details);
}

inline void annotate_tss(Annotations& out, const GeneModel& model, int64_t pos,
                         const AnnotationOptions& options) {
    const auto sites = tss_sites(model, options);
    const auto nearest = nearest_sites(sites, pos, options.max_distance);
    const auto genes = site_genes(nearest, pos, options.max_distance);
    out.strings["TSS"] = labels(genes, options);
    out.strings["TSS_ID"] = ids(genes);
    if (!genes.empty())
        out.integers["TSS_DISTANCE"] = static_cast<int32_t>(magnitude(genes.front().distance));

    if (options.emit.count("neighbors-tss"))
        out.strings["TSS_NEIGHBORS"] = labels(site_genes(sites, pos, options.window, HitOrder::SignedDistance), options, true);
    if (options.emit.count("tss-details")) {
        std::vector<std::string> details;
        for (const auto& site : nearest)
            details.push_back(site_detail(site, options) + (options.canonical_tss ? "|canonical" : "|gene"));
        out.strings["TSS_DETAIL"] = join(details);
    }
    if (options.emit.count("promoter")) annotate_promoters(out, sites, pos, options);
}

inline Annotations annotate(const GeneModel& model, int64_t pos, const AnnotationOptions& options) {
    Annotations out;
    annotate_gene_bodies(out, model, pos, options);
    annotate_neighborhood(out, model, pos, options);
    annotate_tss(out, model, pos, options);
    annotate_features(out, model, pos, options);
    return out;
}
}
