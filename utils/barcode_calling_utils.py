import yaml
import json
import pandas as pd
import math

def queryall(cropped_barcode_dict, query):
    cropped_barcode_list = list(cropped_barcode_dict.keys())

    if query in cropped_barcode_list:
        # is a perfect match
        return 1, cropped_barcode_dict[query]

    else:
        scoredict = {
            sum([1 for x in range(len(query)) if query[x] == y[x]])
            / float(len(query)): y
            for y in cropped_barcode_list
        }
        scores = list(scoredict.keys())
        scores.sort(reverse=True)
        return scores[0], cropped_barcode_dict[scoredict[scores[0]]]
    
def barcodeset(df, barcodecol, genecol):
    barcodeset = {}
    count = 1
    for row in df[[barcodecol, genecol]].itertuples(index=False):
        barcodeset[row[0]] = (count, row[1])
        count += 1
    # output looks like {'GAGTTAGTGAGA': (1, 'SNAPC1')}
    return barcodeset

def match_barcode_to_library(library_location, library_structure, call_col, df):
    library = pd.read_csv(library_location)

    cropped_barcode_dict = {}
    barcodeset_dict = {}
    segment_dict = {}
    count = 0
    for col_to_match in library_structure.keys():
        cropped_barcode_dict[col_to_match] = {
            y[ : library_structure[col_to_match]['n']]: y for y in library[col_to_match]
        }

        segment_dict[col_to_match] = {'start':count,'stop':count+library_structure[col_to_match]['n']}
        count += library_structure[col_to_match]['n']

        barcodeset_dict[col_to_match] = barcodeset(library,col_to_match, library_structure[col_to_match]['namecol'])

    matchdict = {}
    for col_to_match in library_structure.keys():
        scorelist = []
        matchedbarcode = []
        matchedname = []
        matchedid = []
        for eachbarcode in df[call_col]:
            start = segment_dict[col_to_match]['start']
            stop = segment_dict[col_to_match]['stop']
            eachbarcode = eachbarcode[start:stop]

            eachscore, eachmatch = queryall(cropped_barcode_dict[col_to_match], eachbarcode)
            scorelist.append(eachscore)
            matchedbarcode.append(eachmatch)
            matchedname.append(barcodeset_dict[col_to_match][eachmatch][1])
            matchedid.append(barcodeset_dict[col_to_match][eachmatch][0])
        matchdict[col_to_match] = [scorelist, matchedbarcode, matchedname, matchedid]

    return matchdict

def get_spot_quality_dict(method):
    # Ordered best to worst; used to fix category order in plots
    spot_quality_dict = {
        "simple": {
            1: "Perfect",
            2: "Good",
            3: "Acceptable",
            4: "Recombinant",
            5: "BadSameGene",
            6: "BadDiffGene",
        },
        "ignore_short_barcode":{
            1: "Perfect",
            2: "Good",
            3: "Acceptable",
            4: "Bad",
            5: "Recombinant",
            6: "BadLongBar"
        },

    }
    return spot_quality_dict[method]

# assign spot quality categories to each spot if there are multiple score cols
def categorize_spots(df, SBS_score_col, score_cols, library_structure, library_location, spot_quality_method="simple"):
    if spot_quality_method == "simple":
        def do_spot_cats(row):
            score1 = row[score_cols[0]]
            score2 = row[score_cols[1]]

            len1 = library_structure[score_cols[0].replace(f"{SBS_score_col}_",'')]['n']
            len2 = library_structure[score_cols[1].replace(f"{SBS_score_col}_",'')]['n']

            # NOTE - these aren't actual hamming checks, just assuming off-by 1 is correctable
            # We may overcall recom, for example, because in at least some libraries we've handled, the short
            # barcode has only a hamming of 1 and thus is uncorrectable
            hamperfect_score1 = (len1 - 1) / len1
            hamperfect_score2 = (len2 - 1) / len2
            
            gene1 = row[score_cols[0].replace("Score", "GeneCode")]
            gene2 = row[score_cols[1].replace("Score", "GeneCode")]

            if gene1 == gene2:
                if score1 == 1.0 and score2 == 1.0:
                    return "Perfect"
                elif score1 in [1.0, hamperfect_score1] and score2 in [1.0, hamperfect_score2]:
                    return "Good"
                elif (score1 == 1.0 and score2 < hamperfect_score2) or (score1 < hamperfect_score1 and score2 == 1.0):
                    return "Acceptable"
                elif (score1 < hamperfect_score1 or score2 < hamperfect_score2) and (score1 != 1.0 and score2 != 1.0):
                    return "BadSameGene"
                else:
                    return "Uncategorized"
            else:
                if score1 in [1.0, hamperfect_score1] and score2 in [1.0, hamperfect_score2]:
                    return "Recombinant"
                else:
                    return "BadDiffGene"
    elif spot_quality_method == "ignore_short_barcode":
        def determine_max_errors(library_location, library_col_name, len_read):
            def hamming_distance(str1: str, str2: str) -> int:
                """Calculate Hamming distance between two strings."""
                if len(str1) != len(str2):
                    raise ValueError("Strings must be of equal length")
                return sum(c1 != c2 for c1, c2 in zip(str1, str2))

            df = pd.read_csv(library_location)
            full_call = df[library_col_name]
            call_clipper = lambda x: x[:len_read]
            clipped_call = list(full_call.apply(call_clipper))
            hamming_distance_min = len_read
            for x in range(len(clipped_call) - 1):
                for y in range (1, len(clipped_call)):
                    if x != y:
                        hamming_distance_min = min(hamming_distance(clipped_call[x],clipped_call[y]), hamming_distance_min)

            return (math.floor((hamming_distance_min - 1)/2))
            
        def do_spot_cats(row):

            len1 = library_structure[score_cols[0].replace(f"{SBS_score_col}_",'')]['n']
            len2 = library_structure[score_cols[1].replace(f"{SBS_score_col}_",'')]['n']

            long_len = max(len1,len2)
            short_len = min(len1,len2)

            if len1 >len2: #if tie, use first, because hopefully at least less background
                long_index = 0
                short_index = 1
            else:
                long_index = 1
                short_index = 0

            hamming_distance_long = determine_max_errors(library_location, 
                                                            list(library_structure.keys())[long_index],
                                                            long_len)

            hamming_distance_short = determine_max_errors(library_location, 
                                                            list(library_structure.keys())[short_index],
                                                            short_len)

            long_score = row[score_cols][long_index]
            short_score = row[score_cols][short_index]

            min_permissible_score_long = (long_len - hamming_distance_long) /long_len
            min_permissible_score_short = (short_len - hamming_distance_short) /short_len

            gene1 = row[score_cols[0].replace("Score", "GeneCode")]
            gene2 = row[score_cols[1].replace("Score", "GeneCode")]

            if long_score < min_permissible_score_long:
                return "BadLongBar"
            else:
                if gene1 == gene2:
                    if long_score == 1.0 and short_score == 1.0:
                        return "Perfect"
                    elif long_score >= min_permissible_score_long and short_score >= min_permissible_score_short:
                        return "Good"
                    elif long_score >= min_permissible_score_long or short_score < min_permissible_score_short:
                        return "Acceptable"
                else:
                    if long_score >= min_permissible_score_long and short_score >= min_permissible_score_short:
                        return "Recombinant"
                    else:
                        return "Bad"
        # Apply the function. No need for 'args=' if you don't add them to the inner function
    else:
            raise ValueError("Invalid spot quality method specified.")
    
    df['Spot_Category'] = df.apply(do_spot_cats, axis=1)
        
    
        
    assert 'Uncategorized' not in df['Spot_Category'].values, "Error: Some rows were left Uncategorized!"

    return df