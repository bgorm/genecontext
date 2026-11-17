#pragma once
#include "genes.hpp"

namespace gc {
struct Transcript {
    std::string gene_id, id;
    int64_t start = 0, end = 0;
    char strand = '.';
    std::vector<std::string> tags;
    bool tagged(const std::string& tag) const { return std::find(tags.begin(),tags.end(),tag) != tags.end(); }
    std::optional<int64_t> tss() const {
        if (strand == '+') return start;
        if (strand == '-') return end;
        return std::nullopt;
    }
};
struct Segment {
    std::string gene_id, transcript_id;
    FeatureKind kind;
    int64_t start, end;
    auto key() const { return std::tie(gene_id,transcript_id,kind,start,end); }
    bool operator<(const Segment& other) const { return key() < other.key(); }
    bool operator==(const Segment& other) const { return key() == other.key(); }
};
struct GeneModel {
    std::vector<Gene> genes;
    std::vector<Transcript> transcripts;
    std::vector<Segment> segments;
};
// Resolve parents in three passes: eligible genes, their transcripts, then
// exon/CDS links. The indexes live for one query only. A child cannot revive
// a gene excluded by the gene-level filters.
class ModelResolver {
    GeneModel model_;
    std::unordered_map<std::string,std::string> gene_aliases_;
    std::unordered_map<std::string,std::vector<size_t>> transcript_aliases_;
    using TranscriptKey = std::pair<std::string,std::string>; // gene ID, transcript ID
    std::map<TranscriptKey,size_t> transcript_indexes_;

    void add_gene_alias(const std::string& alias, const std::string& id) {
        if (alias.empty()) return;
        const auto [entry, inserted] = gene_aliases_.emplace(alias, id);
        require(inserted || entry->second == id, "Ambiguous annotation identifier: " + alias);
    }

    std::string gene_for(const std::string& id) const {
        const auto found = gene_aliases_.find(id);
        return found == gene_aliases_.end() ? "" : found->second;
    }

    void add_genes(const std::vector<Feature>& features, const GeneFilter& filter,
                   int64_t query_start, int64_t query_end) {
        for (const auto& feature : features) {
            if (feature.kind != FeatureKind::Gene) continue;
            const auto& gene = feature.region;
            if (!filter.accepts(gene) || gene.start > query_end || gene.end < query_start) continue;
            model_.genes.push_back(gene);
            add_gene_alias(gene.id, gene.id);
            add_gene_alias(feature.source_id, gene.id);
        }
        sort_unique(model_.genes);
    }

    void add_transcript(const Feature& feature, const std::string& gene_id) {
        const auto& region = feature.region;
        const TranscriptKey key{gene_id, region.id};
        const auto [entry, inserted] = transcript_indexes_.emplace(key, model_.transcripts.size());
        const size_t index = entry->second;
        if (inserted) {
            model_.transcripts.push_back({gene_id, region.id, region.start, region.end, region.strand, feature.tags});
        } else {
            // Never invent a TSS by choosing an arbitrary transcript part.
            auto& prior = model_.transcripts[index];
            require(prior.start == region.start && prior.end == region.end && prior.strand == region.strand,
                    "Conflicting or multipart transcript coordinates: " + region.id);
            prior.tags.insert(prior.tags.end(), feature.tags.begin(), feature.tags.end());
            sort_unique(prior.tags);
        }
        for (const auto& alias : {region.id, feature.source_id}) {
            if (!alias.empty()) transcript_aliases_[alias].push_back(index);
        }
    }

    void add_transcripts(const std::vector<Feature>& features) {
        for (const auto& feature : features) {
            if (feature.kind != FeatureKind::Transcript) continue;
            std::set<std::string> genes;
            const auto explicit_gene = gene_for(feature.gene_id);
            if (!explicit_gene.empty()) genes.insert(explicit_gene);
            for (const auto& parent : feature.parents) {
                const auto id = gene_for(parent);
                if (!id.empty()) genes.insert(id);
            }
            if (!explicit_gene.empty())
                require(genes.size() == 1, "Conflicting gene_id and Parent for transcript " + feature.region.id);
            for (const auto& gene_id : genes) add_transcript(feature, gene_id);
        }
    }

    void add_segments(const std::vector<Feature>& features) {
        for (const auto& feature : features) {
            if (feature.kind != FeatureKind::Exon && feature.kind != FeatureKind::Cds) continue;
            std::set<TranscriptKey> links;
            const auto explicit_gene = gene_for(feature.gene_id);
            auto link_transcript = [&](const std::string& id) {
                const auto found = transcript_aliases_.find(id);
                if (found == transcript_aliases_.end()) return false;
                for (const size_t index : found->second) {
                    const auto& transcript = model_.transcripts[index];
                    require(explicit_gene.empty() || explicit_gene == transcript.gene_id,
                            "Conflicting gene/transcript relationship at " + feature.region.chrom + ":" +
                            std::to_string(feature.region.start));
                    links.emplace(transcript.gene_id, transcript.id);
                }
                return true;
            };

            // GTF can identify a transcript even when its feature row is absent.
            if (!feature.transcript_id.empty() && !link_transcript(feature.transcript_id) && !explicit_gene.empty())
                links.emplace(explicit_gene, feature.transcript_id);
            for (const auto& parent : feature.parents) {
                if (link_transcript(parent)) continue;
                const auto gene_id = gene_for(parent);
                if (!gene_id.empty()) links.emplace(gene_id, ""); // direct gene-level exon/CDS
            }
            if (links.empty() && !explicit_gene.empty() && feature.parents.empty())
                links.emplace(explicit_gene, "");

            for (const auto& [gene_id, transcript_id] : links)
                model_.segments.push_back({gene_id, transcript_id, feature.kind, feature.region.start, feature.region.end});
        }
        sort_unique(model_.segments);
    }

public:
    static GeneModel resolve(const std::vector<Feature>& features, const GeneFilter& filter,
                             int64_t query_start, int64_t query_end) {
        ModelResolver resolver;
        resolver.add_genes(features, filter, query_start, query_end);
        resolver.add_transcripts(features);
        resolver.add_segments(features);
        return std::move(resolver.model_);
    }
};

inline GeneModel resolve_model(const std::vector<Feature>& features, const GeneFilter& filter,
                               int64_t query_start, int64_t query_end) {
    return ModelResolver::resolve(features, filter, query_start, query_end);
}
}
