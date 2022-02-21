# This script calculates gene similarity using GO and AIC method
# Author:       Joseph Yu
# Date Created: 02-15-2022

import sys
import os
import operator
import mmap
import sample_GO_tree


def cmd_input_validation():

    if len(sys.argv) != 3:
        raise ValueError("Usage: python gene_similarity <gene 1 Ontology Annotation> <gene 2 Ontology Annotation>")

    print(f"Gene 1 Annotation file: {sys.argv[1]}")
    print(f"Gene 2 Annotation file: {sys.argv[2]}")

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    if operator.not_(os.path.exists(file1)) and operator.not_(os.path.exists(file2)):
        raise FileNotFoundError("files do not exist.")

    try:
        file1_content = open(file1)
        file2_content = open(file2)

    except IOError:
        raise(IOError("files cannot be opened"))

    finally:
        file1_content.close()
        file2_content.close()

    return file1, file2


def read_file(file_name):

    information_content = {}
    text = None
    with open(file_name, "r+t") as file_content:
        # memory-map the file, size 0 means whole file
        mm = mmap.mmap(file_content.fileno(), length=0, access=mmap.ACCESS_READ)

        # read from file and place in dictionary
        text = mm.read()


        '''
        # read content via standard file methods
        print(mm.readline())  # prints b"Hello Python!\n"
        # read content via slice notation
        print(mm[:5])  # prints b"Hello"
        # update content using slice notation;
        # note that new content must have same size
        mm[6:] = b" world!\n"
        # ... and read again using standard file methods
        mm.seek(0)
        print(mm.readline())  # prints b"Hello  world!\n"
        '''
        # close the map
        mm.close()

    text = text.decode("utf-8").strip().split(',')

    return information_content


if __name__ == '__main__':
    input1, input2 = cmd_input_validation()

    gene1_information_content = read_file(input1)
    gene2_information_content = read_file(input2)
