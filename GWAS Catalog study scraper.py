import requests
import json
import io
from datetime import datetime

## Converts string into integer if possible
def isInteger(input):
    try:
        int(input)
    except ValueError:
        return False
    else:
        return True

## Returns the manufacturer of each platform joined by commas
def platforms(study):
  manufacturers = [platform["manufacturer"] for platform in study["platforms"]]
  return ','.join(manufacturers)

## Returns discoveryAncestry groups as a comma separated set to remove duplicates
def discoveryAncestry(study):
  discoveryAncestryList = set()
  for ancestry in study["ancestries"]:
    if ancestry["type"] == "initial":
      for group in ancestry["ancestralGroups"]:
        discoveryAncestryList.add(group["ancestralGroup"])
  formattedSet = ', '.join(discoveryAncestryList)
  if formattedSet:
    return formattedSet
  else:
    return "Unknown"
  
## Returns sum of discoverySampleSize
def discoverySampleSize(study):
  ancestries = study["ancestries"]
  discoverySample = sum(ancestry["numberOfIndividuals"] for ancestry in ancestries if ancestry["type"] == "initial" and ancestry["numberOfIndividuals"] is not None)
  return str(discoverySample)

## Returns replicationAncestry groups as a comma separated set to remove duplicates
def replicationAncestry(study):
  replicationAncestryList = set()
  for ancestry in study["ancestries"]:
    if ancestry["type"] == "replication":
      for group in ancestry["ancestralGroups"]:
        replicationAncestryList.add(group["ancestralGroup"])
  formattedSet = ', '.join(replicationAncestryList)
  if formattedSet:
    return formattedSet
  else:
    return "None"
          
## Returns sum of replicationSampleSize
def replicationSampleSize(study):
  ancestries = study["ancestries"]
  replicationSample = sum(ancestry["numberOfIndividuals"] for ancestry in ancestries if ancestry["type"] == "replication" and ancestry["numberOfIndividuals"] is not None)
  return str(replicationSample)

## Returns if study is imputed
def isImputed(study):
  if study["imputed"]:
    return "1"
  else:
    return "0"
  
## Returns a list of chromosomeName counts for both Location Data and rsId location Data
def getChromList(study):
  accessionId = study.get("accessionId")
  if accessionId:
    studyResponse = requests.get(f"https://www.ebi.ac.uk/gwas/rest/api/studies/{accessionId}/snps")
    studyData = studyResponse.json()
    snps = studyData["_embedded"]["singleNucleotidePolymorphisms"]
    # Location Data portion of list
    locationsTally = [0] * 26
    for snp in snps:
      for location in snp["locations"]:
        if isInteger(location["chromosomeName"]):
          locationsTally[int(location["chromosomeName"]) - 1] += 1
        elif location["chromosomeName"] == "X":
          locationsTally[23] += 1
        elif location["chromosomeName"] == "Y":
          locationsTally[24] += 1
        else:
          locationsTally[25] += 1
    # rsId Data portion of list
    rsIdTally = [0] * 26
    for snp in snps:
      if snp["rsId"][0:3] == "chr" and isInteger(snp["rsId"][3:5].split(":")[0]):
          rsIdTally[int(snp["rsId"][3:5].split(":")[0]) - 1] += 1
      elif snp["rsId"][3:5].split(":")[0] == "X":
          rsIdTally[23] += 1
      elif snp["rsId"][3:5].split(":")[0] == "Y":
          rsIdTally[24] += 1
      elif snp["rsId"][0:3] == "chr":
          rsIdTally[25] += 1
    # Merge Location and rsID location lists
    completeList = locationsTally + rsIdTally
    convertToString = [str(number) for number in completeList]
    formattedList = "\t".join(convertToString)
    return formattedList
 
## Creates a .txt file for 1000 study items
def createPage(i):
  ## Create a file to write to with the appropriate name
  fileName = f"Trial_File_Page{i}.txt"
  f = io.open(fileName, "w", encoding="utf-8")
  ## make the API request MAXIMUM size = 108869 *1150*
  api_url = f"https://www.ebi.ac.uk/gwas/rest/api/studies?page={i}&size=1000"
  response = requests.get(api_url)
  ## Parse the JSON into a Python dictionary
  data = response.json()
  ## define studies' location in data
  studies = data["_embedded"]["studies"]
  ## Create document headings
  f.write("accessionId\tpubmedId\ttitle\tpublicationDate\tsnpCount\tplatform\tdiscoveryAncestries\tdiscoverySample\treplicationAncestries\treplicationSample\ttrait\timputed\tsummaryStats\t")
  j = 1
  while j < 24:
    f.write(f"chrom{j}\t")
    j += 1
  f.write("chromX\tchromY\tchromOther\t")
  j = 1
  while j < 24:
    f.write(f"rsIdChrom{j}\t")
    j += 1
  f.write("rsIdChomX\trsIdChomY\trsIdChromOther\n")

  for study in studies:
    ## accessionId
    f.write(study["accessionId"])
    f.write('\t')
    ## pubmedId
    f.write(study["publicationInfo"]["pubmedId"])
    f.write('\t')
    ## studyTitle
    f.write(study["publicationInfo"]["title"])
    f.write('\t')
    ## publicationDate
    f.write(study["publicationInfo"]["publicationDate"])
    f.write('\t')
    ## snpCount
    f.write(str(study["snpCount"]))
    f.write('\t')
    ## platforms
    f.write(platforms(study))
    f.write('\t')
    ## disocveryAncestries
    f.write(discoveryAncestry(study))
    f.write('\t')
    ## discoverySample
    f.write(discoverySampleSize(study))
    f.write('\t')
    ## replicationAncestry
    f.write(replicationAncestry(study))
    f.write('\t')
    ## replicationSample
    f.write(replicationSampleSize(study))
    f.write('\t')
    ## trait
    if study["diseaseTrait"]:
      f.write(study["diseaseTrait"]["trait"])
    else:
      f.write("Unknown")
    f.write('\t')
    ## imputed
    f.write(isImputed(study))
    f.write("\t")
    ## fullPvalueSet
    if study["fullPvalueSet"]:
      f.write("1")
    else:
      f.write("0")
    f.write("\t")
    ## chromosomeList
    f.write(getChromList(study))
    f.write("\n")

  ## Terminate file access
  f.close()
  ## Print when page # complete
  print(f"Page {i} complete!")
  print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

print("Running!")
i = 0
while i < 114:
    createPage(i)
    i += 1