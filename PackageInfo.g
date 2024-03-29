#
# GaloisGroups: Computing Galois Groups using GAP and PARI
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "GaloisGroups",
Subtitle := "Computing Galois Groups using GAP and PARI",
Version := "0.1",
Date := "31/01/2018", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
             rec(
                  IsAuthor := true,
                  IsMaintainer := true,
                  FirstNames := "Bill",
                  LastName := "Allombert",
                  WWWHome := "https://www.math.u-bordeaux.fr/~ballombe/",
                  Email := "Bill.Allombert@math.u-bordeaux.fr",
                  PostalAddress := Concatenation(
                                                 "IMB, UMR XXX",
                                                 "Bâtiment A33"
                                                 ),
                  Place := "Bordeaux",
                  Institution := "CNRS/Université de Bordeaux",
                 ),
             rec(
                  IsAuthor := true,
                  IsMaintainer := true,
                  FirstNames := "Vincent",
                  LastName := "Delecroix",
                  WWWHome := "http://www.labri.fr/perso/vdelecro/",
                  Email := "vincent.delecroix@u-bordeaux.fr",
                  PostalAddress := Concatenation(
                                                 "LaBRI, UMR 5800",
                                                 "Bâtiment A30",
                                                 "351, cours de la Libération",
                                                 "33405 Talence cedex",
                                                 "FRANCE"),
                  Place := "Bordeaux",
                  Institution := "CNRS/Université de Bordeaux",
                 ),
             
             rec(
                  IsAuthor := true,
                  IsMaintainer := true,
                  FirstNames := "Markus",
                  LastName := "Pfeiffer",
                  WWWHome := "https://markusp.morphism.de/",
                  Email := "markus.pfeiffer@st-andrews.ac.uk",
                  PostalAddress := Concatenation(
                                                  "School of Computer Science\n",
                                                  "University of St Andrews\n",
                                                  "Jack Cole Building, North Haugh\n",
                                                  "St Andrews, Fife, KY16 9SX\n",
                                                  "United Kingdom" ),
                  Place := "St Andrews",
                  Institution := "University of St Andrews",
                 ),
],


SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/gap-packages/", ~.PackageName ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
#SupportEmail   := "TODO",
PackageWWWHome  := "https://gap-packages.github.io/GaloisGroups/",
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "GaloisGroups",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Computing Galois Groups using GAP and PARI/GP",
),

Dependencies := rec(
  GAP := ">= 4.11",
  NeededOtherPackages := [ [ "GAPDoc", ">= 1.6" ]
                         , [ "Digraphs", ">= 0.11" ]
                         , [ "ferret", ">= 0.8.0" ]
                         , [ "alnuth", ">= 3.1.0" ]
                         , [ "TransGrp", ">= 2.0.4" ]
                         , [ "PARIInterface", ">= 0.1" ] ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := function()
        return true;
    end,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));


