import yaml
import json
import pandas as pd

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

# assign spot quality categories to each spot if there are multiple score cols
def categorize_spots(df, SBS_score_col, score_cols, library_structure, spot_quality_method="simple"):
    if spot_quality_method == "simple":
        def do_spot_cats(row):
            score1 = row[score_cols[0]]
            score2 = row[score_cols[1]]

            len1 = library_structure[score_cols[0].replace(f"{SBS_score_col}_",'')]['n']
            len2 = library_structure[score_cols[1].replace(f"{SBS_score_col}_",'')]['n']

            hamperfect_score1 = (len1 - 1) / len1
            hamperfect_score2 = (len2 - 1) / len2
            
            gene1 = row[score_cols[0].replace("Score", "GeneCode")]
            gene2 = row[score_cols[1].replace("Score", "GeneCode")]

            if gene1 == gene2:
                if score1 == 1.0 and score2 == 1.0:
                    return "Perfect"
                # 2. Use standard Python 'in' instead of pandas '.isin()' for single floats
                elif score1 in [1.0, hamperfect_score1] and score2 in [1.0, hamperfect_score2]:
                    return "Good"
                elif (score1 == 1.0 and score2 < hamperfect_score2) or (score1 < hamperfect_score1 and score2 == 1.0):
                    return "Acceptable"
                elif (score1 < hamperfect_score1 or score2 < hamperfect_score2) and (score1 != 1.0 and score2 != 1.0):
                    return "Bad"
                else:
                    return "Uncategorized"
            else:
                if score1 in [1.0, hamperfect_score1] and score2 in [1.0, hamperfect_score2]:
                    return "Recombinant"
                else:
                    return "Bad"
                    
        # Apply the function. No need for 'args=' if you don't add them to the inner function
        df['Spot_Category'] = df.apply(do_spot_cats, axis=1)
        
    else:
        raise ValueError("Invalid spot quality method specified.")
        
    assert 'Uncategorized' not in df['Spot_Category'].values, "Error: Some rows were left Uncategorized!"

    return df