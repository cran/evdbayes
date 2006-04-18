".First.lib" <-
function(lib, pkg)
{
  #require(evd)
  library.dynam("evdbayes", package = pkg, lib.loc = lib)
  return(invisible(0))
}

