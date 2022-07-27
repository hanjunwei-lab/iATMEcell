#' @title Get example dataset
#' @description This function is used to achieve example dataset.
#' @param exampleData  A character, should be one of"GEP","clinicaldata",
#' "TMEcellinfo", "TME_related_Goterm", "Jaccardscore" and "GoCellconGene".
#' @return example dataset
#' @usage GetExampleSet(exampleData)
#' @export


GetExampleSet<-function(exampleData){
  if(!exists("envData")) {
    utils::data("envData",package="iATMEcell")
  }

  if (exampleData=="GEP")
  {
    dataset<- get("GEP",envir=envData)
    return(dataset)
  }
  if (exampleData=="clinicaldata")
  {
    dataset<- get("clinicaldata",envir=envData)
    return(dataset)
  }
  if (exampleData=="TMEcellinfo")
  {
    dataset<- get("TMEcellinfo",envir=envData)
    return(dataset)
  }
  if (exampleData=="TME_related_Goterm")
  {
    dataset<- get("TME_related_Goterm",envir=envData)
    return(dataset)
  }
  if (exampleData=="Jaccardscore")
  {
    dataset<- get("Jaccardscore",envir=envData)
    return(dataset)
  }
  if (exampleData=="GoCellconGene")
  {
    dataset<- get("GoCellconGene",envir=envData)
    return(dataset)
  }
  if (exampleData=="iTMEcellresult")
  {
    dataset<- get("iTMEcellresult",envir=envData)
    return(dataset)
  }
  if (exampleData=="GeomSplitViolin")
  {
    dataset<- get("GeomSplitViolin",envir=envData)
    return(dataset)
  }
}
