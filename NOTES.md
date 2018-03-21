Warning Messages (GATK)
================================================================================

### Unk module
    WARN  NativeLibraryLoader - Unable to load libgkl_utils.so from native/libgkl_utils.so (/tmp/nplatt/libgkl_utils6908554598006466719.so: libgomp.so.1: cannot open shared object file: No such file or directory)
    WARN  IntelPairHmm - Intel GKL Utils not loaded

### CombineGVCFs
    WARN  ReferenceConfidenceVariantContextMerger - Detected invalid annotations: When trying to merge variant contexts at location KL250487.1:7703 the annotation MLEAC=[2, 0] was not a numerical value and was ignored

### HaplotypeCaller
    WARN DepthPerSampleHC - Annotation will not be calculated, genotype is not called or alleleLikelihoodMap is null
    WARN PairHMM - ***WARNING: Machine does not have the AVX instruction set support needed for the accelerated AVX PairHmm. Falling back to the MUCH slower LOGLESS_CACHING implementation!
    WARN StrandBiasBySample - Annotation will not be calculated, genotype is not called or alleleLikelihoodMap is null

### GenotypeGVCFs
    WARNING: No valid combination operation found for INFO field MLEAC - the field will NOT be part of INFO fields in the generated VCF records
    WARNING: No valid combination operation found for INFO field MLEAF - the field will NOT be part of INFO fields in the generated VCF records
    WARNING: No valid combination operation found for INFO field DS - the field will NOT be part of INFO fields in the generated VCF records
    WARNING: No valid combination operation found for INFO field DS - the field will NOT be part of INFO fields in the generated VCF records
    WARNING: No valid combination operation found for INFO field InbreedingCoeff - the field will NOT be part of INFO fields in the generated VCF records

### ANY SORTING STEP
htsjdk.samtools.util.RuntimeIOException: java.io.IOException: No space left on device
Caused by: java.io.IOException: No space left on device

These cause the job to fail and seem to be remedied by adding a tmp_dir to the gatk options.  probably something to do with running out of space on the compute nodes themselves


Other
==========
Mitochondira is likely AMPZ01026399.1
