#!/usr/bin/env Rscript

main <- function() {
  if (!dir.exists('assets')) {
    options(timeout=600)
    assets_url <- paste0('https://files.osf.io/v1/resources/ca6zu/providers/osfstorage/',
                         '65cddc3e6c2a400199187f6b/?zip=&_gl=1*l3s9hi*_ga*',
                         'NjQ1MzM5ODAxLjE3MDA2NTIyMTc.*_ga_YE9BMGGWX8*',
                         'MTcwNzk4OTgzNC4yMi4xLjE3MDc5OTAyOTguNjAuMC4w')
    tryCatch({
      cat('Downloading assets from osf: https://osf.io/ca6zu/files/osfstorage')
      download.file(assets_url, 'assets.zip', quiet=T)
      unzip('assets.zip', exdir='assets')
      file.remove('assets.zip')
    }, error=function(e) {
      cat('\n!! ERROR occurred while downloading assets from OSF:\n!!', conditionMessage(e))
      cat('\n!! Try manual download from:\n!!', assets_url)
      return(invisible())
    })
  }
  cat('\nAll required assets downloaded and available!')
  cat('\n', rep('=', 45), sep='')
}

main()