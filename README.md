# AETG
Simple two implementation of AETG algorithm

Here i give two different implementation with similar ideas.Both implementation use C++ language. With two different purpose,
the former purpose is to generate minimum testcase as less time as possible, the latter purpose is also to generate minimum 
testcase without putting time savings first.

We call the first AETG1 and the second AETG2 and both source code can be found in branch src,which supports t-way
combinatorial test generation for up to 10-way coverage as i take usual circumstance into consideration. But you can achieve n-way coverage by modifying the source code a little bit.

Notice that in addition to generating minimum test case, i dont achieve another function which is frequently talked in combinatorial testing,such as constriant, seed and mixed strength, ect.

