# This script calculates gene similarity using GO and AIC method
# Author:       Joseph Yu
# Date Created: 02-15-2022

import sys
import os
import operator
import mmap
import math

import sample_GO_tree
import annotation_data

information_content_dict = {}


def cmd_input_validation():
    if len(sys.argv) != 3:
        raise ValueError("Usage: python gene_similarity <gene 1 Ontology Annotation> <gene 2 Ontology Annotation>")

    print(f"Gene 1 Annotation file: {sys.argv[1]}")
    print(f"Gene 2 Annotation file: {sys.argv[2]}")

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    if operator.not_(os.path.exists(file1)) and operator.not_(os.path.exists(file2)):
        raise FileNotFoundError("files do not exist.")

    file1_content = None
    file2_content = None
    try:
        file1_content = open(file1)
        file2_content = open(file2)

    except IOError:
        raise(IOError("files cannot be opened"))

    finally:
        file1_content.close()
        file2_content.close()

    return file1, file2


# @TODO return actual annotation data
def get_annotation(term):
    return annotation_data.annotation_data_dict[term][0]


def get_annotation_type(term):
    return annotation_data.annotation_data_dict[term][1]


def traverse(term, tree, frequency_result, finished):
    if tree[term][0] == "None":
        if term not in frequency_result.keys():
            frequency_result[term] = int(get_annotation(term))
            finished.add(term)
    else:
        for child in tree[term]:
            if term not in frequency_result.keys():
                frequency_result[term] = int(get_annotation(term))
            # print(f"{term} {child} {frequency_result}")
            if child not in finished:
                traverse(child, tree, frequency_result, finished)
                finished.add(child)
            frequency_result[term] += frequency_result[child]


# @TODO function for calculating the IC from data
def calculate_information_content(tree, frequency_result):
    # @TODO traverse tree and calculate IC
    finished = set()
    traverse("GO:0008150", tree, frequency_result, finished)
    traverse("GO:0003674", tree, frequency_result, finished)

    bp_root = frequency_result["GO:0008150"]
    ml_root = frequency_result["GO:0003674"]
    for term in frequency_result:
        try:
            if get_annotation_type(term) == "BP":
                frequency_result[term] = (-1) * math.log(frequency_result[term] / bp_root, 10)
            else:
                frequency_result[term] = (-1) * math.log(frequency_result[term] / ml_root, 10)
        except ZeroDivisionError:
            if get_annotation_type(term) == "BP":
                print(frequency_result[term], " ", frequency_result["GO:0008150"])
            else:
                print(frequency_result[term], " ", frequency_result["GO:0003674"])


def get_information_content(term):
    return information_content_dict[term]


def get_ancestor(term, tree, ancestors):
    ancestors.add(term)
    for ancestor in tree[term]:
        if ancestor == "None":
            return
        else:
            get_ancestor(ancestor, tree, ancestors)


def read_file(file_name):
    text = None
    with open(file_name, "r+t") as file_content:
        # memory-map the file, size 0 means whole file
        mm = mmap.mmap(file_content.fileno(), length=0, access=mmap.ACCESS_READ)

        # read from file and place in dictionary
        text = mm.read()

        # close the map
        mm.close()

    text = text.decode("utf-8").strip().split(',')
    print(f"{file_name} Ontology Annotation: {text}")

    return text


def list_to_information_content(terms):
    information_content = {}
    # Adding Information Content Values
    for term in terms:
        information_content[term] = get_information_content(term)

    return information_content


def to_knowledge(information_content):
    return 1 / information_content


def to_semantic_weight(knowledge):
    return float(1 / (1 + pow(math.e, (-1) * knowledge)))


def list_to_knowledge(information_content):
    for key in information_content:
        value = information_content[key]
        information_content[key] = to_knowledge(value)
    return information_content


def list_to_semantic_weight(knowledge):
    for key in knowledge:
        value = knowledge[key]
        knowledge[key] = to_semantic_weight(value)
    return knowledge


def gene_sim(ontology_m, m, n):
    polynomial = 1 / (m * n)
    summation = 0
    for row in ontology_m:
        for col in row:
            summation += col
    return polynomial * summation


def main():
    input1, input2 = cmd_input_validation()
    print()  # Print new line separation

    global information_content_dict
    calculate_information_content(sample_GO_tree.GO_tree_sample, information_content_dict)
    terms1 = read_file(input1)
    terms2 = read_file(input2)

    # find ancestor terms of ontology annotations
    ancestor_dict = {}

    # find sets of ancestors for each annotation
    # @TODO remove sample data
    for k in terms1:
        ancestor_set = set()
        get_ancestor(k, sample_GO_tree.GO_tree_sample_reversed, ancestor_set)
        ancestor_dict[k] = ancestor_set
    for k in terms2:
        ancestor_set = set()
        if k not in ancestor_dict.keys():
            get_ancestor(k, sample_GO_tree.GO_tree_sample_reversed, ancestor_set)
            ancestor_dict[k] = ancestor_set

    gene1_information_content = list_to_information_content(terms1)
    gene2_information_content = list_to_information_content(terms2)

    # convert gene annotation ontology to knowledge
    gene1_knowledge = list_to_knowledge(gene1_information_content)
    gene2_knowledge = list_to_knowledge(gene2_information_content)

    # convert gene annotation ontology knowledge to semantic weight
    gene1_sw = list_to_semantic_weight(gene1_knowledge)
    gene2_sw = list_to_semantic_weight(gene2_knowledge)

    # find semantic weight from all ancestor terms
    sw_dictionary = {}
    for k in ancestor_dict:
        for element in ancestor_dict[k]:
            if element == "GO:0008150" or element == "GO:0003674":
                sw_dictionary[element] = 0
            else:
                sw_dictionary[element] = to_semantic_weight(to_knowledge(get_information_content(element)))

    # calculate semantic value for gene1 and gene2 ontology annotations
    # from semantic weight of all ancestors (inclusive)
    gene1_sv = {}
    gene2_sv = {}
    for gene in gene1_sw:
        sum_sv = 0
        term_set = ancestor_dict[gene]
        for term in term_set:
            sum_sv += sw_dictionary[term]

        gene1_sv[gene] = sum_sv

    for gene in gene2_sw:
        sum_sv = 0
        term_set = ancestor_dict[gene]
        for term in term_set:
            sum_sv += sw_dictionary[term]

        gene2_sv[gene] = sum_sv

    # calculate similarity between gene1 and gene2 ontology annotations
    ontology_matrix = [[None for x in range(len(gene2_sv))] for y in range(len(gene1_sv))]

    gene1_counter = 0
    for gene1_ont in gene1_sv:
        gene2_counter = 0
        for gene2_ont in gene2_sv:
            sum_semantic_weights = 0
            # add up semantic weights of all common ancestors
            common_ancestors = ancestor_dict[gene1_ont].intersection(ancestor_dict[gene2_ont])
            for ancestor_term in common_ancestors:
                sum_semantic_weights += sw_dictionary[ancestor_term]

            sim = sum_semantic_weights / (gene1_sv[gene1_ont] + gene2_sw[gene2_ont])
            ontology_matrix[gene1_counter][gene2_counter] = sim
            gene2_counter += 1
        gene1_counter += 1

    # calculate gene similarities
    score = gene_sim(ontology_matrix, len(gene1_sv), len(gene2_sv))
    print(f"Gene similarity score is: {score}")


if __name__ == '__main__':
    main()
