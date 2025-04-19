print("Testing Snakemake interface")
print(exists("snakemake"))
if(exists("snakemake")) {
  print(names(snakemake))
}
print("Done testing")