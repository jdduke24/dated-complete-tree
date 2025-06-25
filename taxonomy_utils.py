import re

tx_levels = {
                # need some value for these, though they should be overwritten during loading and annotation depending on
                # whether in plants (plant[sub]section) or animals (animal[sub]section)
                'section': 999,
                'subsection': 999,
                # no information, could be anywhere, assume not superior to anything
                'no rank':-2,
                'mrca':-1,
                # these are leaf ranks - all below species
                'no rank - terminal':1,
                'varietas':2,
                'variety':2,
                'forma':2,
                'infraspecificname':2,
                'subspecies':2,
                # species
                'species':3, 'species (promoted)':3, 'species (imputed)':3, 'species subgroup':4, 'species group':5,
                # section in plants
                'plantsubsection':6, 'plantsection':7,
                # genus
                'subgenus':8, 'genus':9,
                # tribe
                'subtribe':10, 'tribe':11, 'supertribe':12,
                # family
                'subfamily':13, 'family':14, 'superfamily':15,
                # section in animals
                'animalsubsection':16, 'animalsection':17,
                # order
                'parvorder':18, 'infraorder':19, 'suborder':20, 'order':21, 'superorder':22,
                # cohort
                'subcohort':23, 'cohort':24,
                # class
                'infraclass':25, 'subclass':26, 'class':27, 'superclass':28,
                # phylum
                'infraphylum':29, 'subphylum':30, 'phylum':31, 'superphylum':32,
                # kingdom
                'infrakingdom':33, 'subkingdom':34, 'kingdom':35,
                # domain
                'domain':36
            }


name_regex = []

normal_name = r"[A-Za-z×ëäöﬂæ«\-\[\]]"

anything = r"[A-Za-z0-9αβγδεμθ\.,/&#?=:()<>+*_\[\]\-]"

ottid = r"(ott[0-9]+)"

# 0. standard name with x (extinct): x_Genus_species_ott1234. Mark as extinct and ignore.
name_regex.append("^x_(" + normal_name + "+)_(" + normal_name + "+)_" + ottid + "$")

# 1. uncultured with sp. at end: uncultured_BL#AH_sp._ott1234. Use whole name as species; whole name without sp. as genus.
name_regex.append("^[Uu]ncultured_(?:Candidatus_)?(" + anything + r"+)_(sp\.)_" + ottid + "$")

# 2. uncultured without sp. at end: uncultured_BL#AH_ott1234. Use whole name as both species and genus.
name_regex.append("^[Uu]ncultured_(" + anything + "+)_" + ottid + "$")

# 3. unidentified: unidentified_BL#AH_ott1234. Use whole name as both species and genus.
name_regex.append("^[Uu]nidentified_(" + anything + "+)_" + ottid + "$")

# 4. standard name with Candidatus: Candidatus_Genus_species_ott1234. Use genus name and species name as given.
name_regex.append("^[Cc]andidatus_(" + normal_name + "+)_(" + normal_name + "+)_" + ottid + "$")

# 5. weird name with Candidatus: Candidatus_BLAH#AS_121_ott1234. Use whole name as both species and genus.
name_regex.append("^[Cc]andidatus_(" + anything + "+)_" + ottid + "$")

# 6. anything with cf., sp., aff., nr. at start of name. Use normal name part as genus, whole name without cf as species.
name_regex.append(r"^(cf\.|sp\.|aff\.|nr\.)_(" + normal_name + "+)_(" + anything + "+)_" + ottid + "$")

# 7. any other name that starts with lower case: use whole name as genus and species ("do no harm")
name_regex.append("^([a-z].*)_" + ottid + "$")

# 8. genus only with sp. at end: Genus_sp._ott1234. Use whole name as species; whole name without sp. as genus.
name_regex.append("^(" + normal_name + r"+)_(sp\.)_" + ottid + "$")

# 9. standard name: Genus_species_ott1234. Use genus name and species name as given.
name_regex.append("^(" + normal_name + "+)_(" + normal_name + "+)_" + ottid + "$")

# 10. with subspecies: Genus_species_subspecies_ott1234. Use genus name and species name as given.
name_regex.append("^(" + normal_name + "+)_(" + normal_name + "+)_(" + normal_name + "+)_" + ottid + "$")

# 11. with variety/forma/subsp: Genus_species_var._subspecies_ott1234. Use genus name and species name as given.
name_regex.append("^(" + normal_name + "+)_(" + normal_name + r"+)_(var\.|f\.|subsp\.)+_(" + anything + "+)_" + ottid + "$")

# 12. genus only: Genus_ott1234. Use genus name as both species and genus
name_regex.append("^(" + normal_name + "+)_" + ottid + "$")

# 13. anything with cf., sp., aff., nr. in name. Use whole name as species; whole name without sp. and any subsequent stuff as genus.
name_regex.append("^(" + anything + r"+)_(cf\.|sp\.|aff\.|nr\.|sect\.|subsect\.|subgen\.|ser\.|trib\.|str\.)(" + anything + "+)_" + ottid + "$")

# 14. with some non-basic-subspecies text: Genus_species_2VRR_BL#AH_ott1234
name_regex.append("^(" + normal_name + "+)_(" + normal_name + "+)_(" + anything + "+)_" + ottid + "$")

# 15. with some non-basic-species text: Genus_2VRR_BL#AH_ott1234
name_regex.append("^(" + normal_name + "+)_(" + anything + "+)_" + ottid + "$")

# 16. with some non-basic-Genus text: BL#AH_ott1234
name_regex.append("^(" + anything + "+)_" + ottid + "$")

standard_name = set([1,4,8,9,10,11,14,15])

single_name = set([2,3,5,7,12,16])

cf_name = set([13])

initial_cf_name = set([6])

def get_genus_and_species(name, ignore_extinct=True):
    """Take a node name for a species (or species subgroup or group) from the Open Tree of Life, and
    return a tuple with our best guess for the genus name and species name.
    """

    for i, regex in enumerate(name_regex):
        m = re.search(name_regex[i], name)
        if m:
            if i == 0 and ignore_extinct:
                genus = "extinct"
                species = "extinct"
                break
            elif i == 0: # include extinct
                genus = m.group(1)
                species = m.group(1) + " " + m.group(2)
                break
            elif i in standard_name:
                genus = m.group(1)
                species = m.group(1) + " " + m.group(2)
                break
            elif i in single_name:
                genus =  m.group(1)
                species = m.group(1)
                break
            elif i in cf_name:
                genus =  m.group(1)
                species = m.group(1) + " " + m.group(2) + m.group(3)
                break
            elif i in initial_cf_name:
                genus =  m.group(2)
                species = m.group(2) + " " + m.group(3)
                break

    return (genus, species)
