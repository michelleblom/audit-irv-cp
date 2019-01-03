# Reference implementation code for pseudo-random sampler
# for election audits or other purposes.
# Written by Ronald L. Rivest
# filename: sampler.py
# url: http://people.csail.mit.edu/rivest/sampler.py
sampler_version = "November 14, 2011"
# 
# Relevant to document being produced by an ad-hoc working group chaired
# by Prof. Philip Stark (U.C. Berkeley) regarding election auditing.
# Tested using python version 2.6.7.   (see www.python.org)
# (Will not work with Python version 3, e.g. 3.x.y)
# (Note added 2014-09-07: As per a suggestion by Chris Jerdonek, one should
#  consider this proposal as based on the use of  UTF-8 encoding for strings 
#  throughout.  This comment resolves some potential ambiguities about how 
#  strings are converted to byte sequences before hashing, and the types of
#  strings input by raw_input, etc.  See
#    https://github.com/cjerdonek/rivest-sampler-tests
# for more discussion and test-cases.
# )


# Edited by Michelle Blom for use in auditing IRV elections
# Feb 2018

# $ python sampler.py seed #ballots

################################################################################
## Standard "MIT License"  http://www.opensource.org/licenses/mit-license.php ##
################################################################################
"""
Copyright (c) 2011 Ronald L. Rivest

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
################################################################################
################################################################################

# import library of cryptography hash functions
# This program uses SHA-256 hash function
# For reference, see, e.g.
#     http://en.wikipedia.org/wiki/SHA-2
#     http://csrc.nist.gov/publications/fips/fips180-2/fips180-2.pdf
# SHA-256 implemented here as hashlib.sha256
import hashlib                       

# import library of string-related functions
import string

import sys

def generate_outputs(n,a,b,seed):
    # check that input parameters are valid
    assert n >= 0
    assert a <= b
    N = (b - a + 1)                        # size of set to draw from 

    #initialization
    output_list = [ ]
    count = 0        

    # loop until we have generated the desired sample of size n
    while count < n:
        count = count + 1
        
        # hash_input is seed followed by comma followed by decimal rep of count
        hash_input = bytes(seed + "," + str(count), "utf-8")

        # Apply SHA-256, interpreting hex output as hexadecimal integer
        # to yield 256-bit integer (a python "long integer")
        hash_output = int(hashlib.sha256(hash_input).hexdigest(),16)
        # determine "pick" as pseudo-random value in range a to b, inclusive,
        # as a function of hash_output
        pick = int(a + (hash_output % (b-a+1)))
       
        output_list.append(pick) 

    return output_list

seed = sys.argv[1]
n = int(sys.argv[2])

plist = generate_outputs(n, 0, n-1, seed)
for i in range(n):
    print(plist[i],end='')
    if i != n-1:
        print(",", end='')

print("\n", end='')
