KL_single <- function(nreps,
                      nsites,
                      mu,
                      phi,
                      theta,
                      sfx,
                      sfy,
                      sd_multi,
                      seed,
                      Candidates) {
    res.sim <- SimDat(
        nreps = nreps,
        nsites = nsites,
        mu = mu,
        phi = phi,
        theta = theta,
        sfx = sfx,
        sfy = sfy,
        sd_multi,
        seed = seed
    )

    temp_kl <- CompareSimRealCount(
        Sim.Candidates = res.sim,
        Real.Candidates = Candidates,
        PLOT = FALSE
    )
    return(temp_kl)
}
