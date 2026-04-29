using CSV
using DataFrames
using FASTX

fasta_file = "data/Genes/gg_13_5_ssualign.fasta"
taxonomy_file = "data/Genes/gg_13_5_taxonomy.txt"
output_seqs = "data/Genes/gg_13_5_ssualign_filtered_30000.txt"
output_labels = "data/Genes/gg_13_5_labels_all_30000.csv"

max_sequences = 30000

num_levels = 7
level_names = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

level_dicts = [Dict{String, Int}() for _ in 1:num_levels] 

level_counters = ones(Int, num_levels)

tax_df = CSV.read(taxonomy_file, DataFrame, delim='\t', header=false, types=String)

tax_map = Dict(row[1] => row[2] for row in eachrow(tax_df))

println("Loaded $(length(tax_map)) taxonomy definitions.")
println("Starting extraction of up to $max_sequences sequences...")
open(output_seqs, "w") do out_seq_f
    open(output_labels, "w") do out_label_f
        
        println(out_label_f, join(level_names, ","))
        
        reader = FASTA.Reader(open(fasta_file, "r"))
        
        seq_count = 0 
        current_labels = fill(-1, num_levels) # placeholder for the 7  labels
        
        for record in reader
            if seq_count >= max_sequences
                break
            end
            
            # Get ID and the actual sequence
            seq_id = FASTA.identifier(record)
            seq = FASTA.sequence(String, record)
            
            # look taxonomy text 
            full_taxonomy = get(tax_map, seq_id, "unknown")
            
            if full_taxonomy == "unknown"
                println("Warning: ID $seq_id not found in taxonomy map.")
                tax_levels = ["unknown"]
            else
                # Split taxonomy test into categories
                tax_levels = split(full_taxonomy, "; ")
            end
            
            # Transform text categories to integer labels for all 7 levels
            for i in 1:num_levels
                limit = min(i, length(tax_levels))
                
                path = join(tax_levels[1:limit], "; ")
                
                if !haskey(level_dicts[i], path)
                    level_dicts[i][path] = level_counters[i]
                    level_counters[i] += 1
                end
                
                current_labels[i] = level_dicts[i][path]
            end
            
            println(out_seq_f, seq)
            
            println(out_label_f, join(current_labels, ","))
            
            seq_count += 1
        end
        
        close(reader)
        
        println("\nFinished extracting $seq_count sequences.")
        
    end
end

println("Created $(output_seqs)")
println("Created $(output_labels)")
println("\nTotal unique clusters found per level:")

for i in 1:num_levels
    println(rpad(level_names[i] * ":", 10), level_counters[i] - 1)
end