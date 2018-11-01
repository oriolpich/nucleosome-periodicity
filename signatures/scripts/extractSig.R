
args = commandArgs(trailingOnly=TRUE)

library('sigfit')

signature_fitting <-function(input, output_data){

    data_read <-read.table(input, header = FALSE,sep="\t")

    # add header and remove unwanted columns
    colnames(data_read) <- c('Sample', 'Chr', 'Pos', 'Ref', 'Alt', 'Trinuc' )
    keeps <- c('Sample','Ref', 'Alt', 'Trinuc')
    data_read <- data_read[keeps]

    # feed it into sigfit
    counts_data <- build_catalogues(data_read)


    data("cosmic_signatures", package = "sigfit")

    # run signatures fitting
    mcmc_samples_fit <- sigfit::fit_signatures(counts = counts_data,
                                               signatures = cosmic_signatures,
                                               iter = 13000,
                                               warmup = 3000,
                                               chains = 1,
                                               seed = 1)
    # get exposures
    exposures <- retrieve_pars(mcmc_samples_fit,
                             feature = "exposures",
                             hpd_prob = 0.90)

    # add COSMIC names
    colnames(exposures$mean) <- rownames(cosmic_signatures)

    # add sample names
    n = data.matrix(unique(data_read$Sample) )
    colnames(n) <-'Sample'

    all_m = cbind(exposures$mean, n)

    write.table(all_m, output_data, sep ="\t",
                row.names = FALSE)

}

signature_fitting(args[1], args[2])