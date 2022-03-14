GO_tree_sample_reversed = {"GO:0022904": ["GO:0045333", "GO:0022900"],
                           "GO:0045333": ["GO:0015980"],
                           "GO:0015980": ["GO:0006091"],
                           "GO:0022900": ["GO:0006091"],
                           "GO:0006091": ["GO:0044237"],
                           "GO:0044237": ["GO:0009987", "GO:0008152"],
                           "GO:0009987": ["GO:0008150"],
                           "GO:0008152": ["GO:0008150"],
                           "GO:0008270": ["GO:0046914"],
                           "GO:0046914": ["GO:0046872"],
                           "GO:0046872": ["GO:0043169"],
                           "GO:0043169": ["GO:0043167"],
                           "GO:0043167": ["GO:0005488"],
                           "GO:0005488": ["GO:0003674"],
                           "GO:0020037": ["GO:0046906"],
                           "GO:0046906": ["GO:1901363", "GO:0097159"],
                           "GO:1901363": ["GO:0005488"],
                           "GO:0097159": ["GO:0005488"],
                           "GO:0005515": ["GO:0005488"],
                           "GO:0004023": ["GO:0004022"],
                           "GO:0004024": ["GO:0004022"],
                           "GO:0004022": ["GO:0018455"],
                           "GO:0018455": ["GO:0016616"],
                           "GO:0016616": ["GO:0016614"],
                           "GO:0016614": ["GO:0016491"],
                           "GO:0016491": ["GO:0003824"],
                           "GO:0003824": ["GO:0003674"],
                           "GO:0004174": ["GO:0022904", "GO:0009055", "GO:0016649"],
                           "GO:0009055": ["GO:0022900", "GO:0016645"],
                           "GO:0016649": ["GO:0016645"],
                           "GO:0016645": ["GO:0016491"],
                           "GO:0008150": ["None"],
                           "GO:0003674": ["None"]
                           }

# Reverse Tree Edges
GO_tree_sample = {}
for value in GO_tree_sample_reversed.keys():
    for term in GO_tree_sample_reversed[value]:
        if term in GO_tree_sample.keys():
            GO_tree_sample[term].append(value)
        else:
            GO_tree_sample[term] = [value]

GO_tree_sample["GO:0004174"] = ["None"]
GO_tree_sample["GO:0004023"] = ["None"]
GO_tree_sample["GO:0004024"] = ["None"]
GO_tree_sample["GO:0005515"] = ["None"]
GO_tree_sample["GO:0020037"] = ["None"]
GO_tree_sample["GO:0008270"] = ["None"]
