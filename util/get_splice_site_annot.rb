#!/usr/bin/env ruby

# get all splice sites (intron/exon boundaries) from UCSC knownGene table

if ARGV.size < 1
  $stderr.puts "USAGE: #{$0} knownGene.txt"
  exit 1
end

kg_fn = ARGV.shift

File.open(kg_fn).each_line do |line|
  tx_id, chr, strand, a, b, c, d, n_exons, exon_start_list, exon_stop_list = 
    line.chomp.split(/\t/)

  n_exons = n_exons.to_i
  if n_exons > 1
    exon_starts = exon_start_list.split(/,/).map { |s| s.to_i }
    exon_stops = exon_stop_list.split(/,/).map { |s| s.to_i }
    # for each internal exon
    1.upto(n_exons-2) do |i|
      puts "#{chr}\t#{exon_starts[i]}\t#{1+exon_starts[i]}\tss\t0\t#{strand}"
      puts "#{chr}\t#{exon_stops[i]}\t#{1+exon_stops[i]}\tss\t0\t#{strand}"
    end
  end
end
