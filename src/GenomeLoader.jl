module GenomeLoader

using BioSequences, GenomicFeatures, FASTX

export load_genome, get_seq


struct GenomeRecord{T}
	name::String
	metadata::String
	seq::T
end

function convert_N_to_gap!(seq)
    for i = 1:length(seq)
        (seq[i] == DNA_N) && (seq[i] = DNA_Gap)
    end
end

function load_genome(fastafile; verbose=true)
	indexfile = fastafile*".fai"
	!isfile(indexfile) && error("Index: $indexfile not found")
	reader = FASTA.Reader(open(fastafile, "r"), index=indexfile)

    genome = Dict{String, GenomeRecord{LongSequence{DNAAlphabet{4}}}}()

    for record in reader
        id = FASTA.identifier(record)
        verbose && println("Loading $id...")
        seq = FASTA.sequence(record)
        convert_N_to_gap!(seq)
        genome[id] = GenomeRecord(id, "", seq)
    end
	genome
end


@inline function get_seq(genome, chr, range)
    ((range.start < 1) || (range.stop .> length(genome[chr].seq))) && return genome[chr].seq[1:0]
    genome[chr].seq[range]
end

@inline function get_seq(genome, chr, range, strand::String)
    seq = get_seq(genome, chr, range)
    (strand == "-") && reverse_complement!(seq)
    seq
end

@inline function get_seq(genome, chr, range, strand::Strand)
    seq = get_seq(genome, chr, range)
    (strand == STRAND_NEG) && reverse_complement!(seq)
    seq
end



end # module
