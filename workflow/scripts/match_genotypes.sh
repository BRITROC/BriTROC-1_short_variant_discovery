# A script which filters VCF files so that all samples have matching predicted genotypes
# The purpose of this is to remove spurious variant calls in which the called variant is not present in all technical replicates

vcf_file=$1

java -jar ../jvarkit/dist/vcffilterjs.jar -e 'function accept(vc) {for(var i=1;i<vc.getNSamples();i++) if(!vc.getGenotype(0).sameGenotype(vc.getGenotype(i))) return false; return true;}accept(variant); ' $vcf_file
