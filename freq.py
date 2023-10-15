import hail as hl
import matplotlib.pyplot as plt

# Initialize Hail
hl.init()

# Load the VCF file
vcf_file = "gnomad.vcf"
mt = hl.import_vcf(vcf_file)

# Calculate allele frequencies
mt = mt.annotate_rows(AF=hl.agg.fraction(hl.agg.collect(mt.GT.n_alt_alleles())))

# Collect the top 10 variants with the highest allele frequencies
top_variants = mt.rows().select("AF", "locus", "alleles").order_by(hl.desc("AF")).head(10)

# Print the top variants
print("Top 10 Variants with Highest Allele Frequencies:")
for i, variant in enumerate(top_variants):
    print(f"{i+1}. {variant.alleles} - Allele Frequency: {variant.AF:.4f}")

# Plot allele frequencies
allele_frequencies = mt.rows().select("AF").collect()
plt.hist(allele_frequencies, bins=20, edgecolor='black')
plt.xlabel("Allele Frequency")
plt.ylabel("Count")
plt.title("Distribution of Allele Frequencies")
plt.show()

# Stop Hail
hl.stop()
