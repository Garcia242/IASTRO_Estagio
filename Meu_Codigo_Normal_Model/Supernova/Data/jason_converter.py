# code to convert txt files into json
# This seems to be quite useful https://www.geeksforgeeks.org/convert-text-file-to-json-in-python/

import json

filename = "supernova.txt"

dict1 = {}

with open(filename) as fh:
    
    for line in fh:
        
        command, description = row.strip().split(None, 1)
        
        dict1[command] = description.strip()
        
out_file = open("supernova1.json", "w")
json.dump(dict1, out_file, indent = 4, sort_keys = False)
out_file.close()
