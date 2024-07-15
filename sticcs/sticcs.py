#!/usr/bin/env python

import numpy as np
from collections import defaultdict
import itertools
import argparse
import sys, gzip, string
import cyvcf2
from sticcs.dac import *

def get_patterns_and_matches(der_counts, exclude_singletons=False):
    patterns, matches, n_matches = np.unique(der_counts, axis=0, return_inverse=True, return_counts=True)
    
    #reorder by abundance
    order = np.argsort(n_matches)[::-1]
    
    n_matches = n_matches[order]
    
    patterns=patterns[order]
    
    conversion_dict = dict(zip(order,range(len(order))))
    
    #efficiently translate the matches to the new order
    #got this from https://stackoverflow.com/questions/16992713/translate-every-element-in-numpy-array-according-to-key
    matches = np.vectorize(conversion_dict.__getitem__)(matches)
    
    if exclude_singletons:
        cutoff = np.where(n_matches > 1)[0][-1] 
        patterns = patterns[0:(cutoff+1)]
        n_matches = n_matches[0:(cutoff+1)]
        matches = matches[matches <= cutoff]
    
    return (patterns,matches,n_matches)


def passes_derived_gamete_test(pattern0, pattern1, ploidies):
    #n_ind = len(ploidies)
    #assert len(pattern0) == len(pattern1), "Number of counts must be the same for both patterns"
    #assert len(pattern0) == n_ind, "Number of counts must match number of ploidy values"
    
    _diff_ = pattern0 - pattern1
    _sum_ = pattern0 + pattern1
    
    #if a position has more derived alleles than the other in an individual, then we know that the 1--0 haplotype exists
    #If we see the inverse in another individual, then we know that the 0--1 haploype also exists
    #if, in addition, there is an individual with more derived alleles at both sites than can be explained
    # by only 0--1 and 1--0 haploypes, then we know that the 1--1 haplotype exists.
    if np.any(_diff_ < 0) & np.any(_diff_ > 0) & np.any(_sum_ > ploidies):
        return False
    
    return True


def compatible_with_all(queryID, otherIDs, patterns, ploidies):
    for otherID in otherIDs:
        if not passes_derived_gamete_test(patterns[queryID], patterns[otherID], ploidies): return False
    
    return True


def get_clusters(patterns, matches, positions, ploidies, second_chances=False, seq_start=None, seq_len=None, silent=False):
    
    total=len(matches)
    
    if not silent: print("\n    Making base clusters.\n    ", file=sys.stderr)
    
    onePercent = int(total/100)
    
    # Make the simplesty possible clusters, with one SNP in each
    # We make two versions of this, because we'll later import SNPs from left and right separately
    clusters_forward =  [{"IDs":set([matches[i]]), "nearest":{matches[i]:positions[i]}, "n":{matches[i]:1}} for i in range(total)]
    clusters_reverse =  [{"IDs":set([matches[i]]), "nearest":{matches[i]:positions[i]}, "n":{matches[i]:1}} for i in range(total)]
    
    if second_chances: knockouts = {"IDs":set(), "nearest":{}, "n":{}, "by":{}}
    
    if not silent: print("\n    Forward search for compatible SNPs\n    ", file=sys.stderr)
    
    #now try importing SNPs from let and right (if compatible) and tracking their nearest position
    for i in range(1,total):
        if not silent and i % onePercent == 0: print(".", end="", file=sys.stderr, flush=True)
        
        #if this SNP is in the previous cluster, import everything
        if matches[i] in clusters_forward[i-1]["IDs"]:
            clusters_forward[i]["IDs"].update(clusters_forward[i-1]["IDs"])
            clusters_forward[i]["nearest"].update(clusters_forward[i-1]["nearest"])
            clusters_forward[i]["nearest"][matches[i]] = positions[i]
            clusters_forward[i]["n"].update(clusters_forward[i-1]["n"])
            clusters_forward[i]["n"][matches[i]] += 1
        
        #if current SNP pattern is not in the previous cluster, only import compatibles and send knockouts to the knockout dictionary
        else:
            if second_chances: new_knockout_IDs = set()
            for ID in clusters_forward[i-1]["IDs"]:
                if passes_derived_gamete_test(patterns[ID], patterns[matches[i]], ploidies):
                    clusters_forward[i]["IDs"].add(ID)
                    clusters_forward[i]["nearest"][ID] = clusters_forward[i-1]["nearest"][ID]
                    clusters_forward[i]["n"][ID] = clusters_forward[i-1]["n"][ID]
                elif second_chances:
                    new_knockout_IDs.add(ID)
            
            if second_chances:
                #check the knockouts in case some can be brought back
                second_chance_IDs = set()
                no_more_chances_IDs = set()
                for ID in knockouts["IDs"]:
                    #first see if it's compatible with the new SNP
                    if passes_derived_gamete_test(patterns[ID], patterns[matches[i]], ploidies):
                        #now see if its nemesis has been knocked out
                        if knockouts["by"][ID] in new_knockout_IDs:
                            #if we get here, this knockout should be compatible with the new cluster, so give it a second chance
                            second_chance_IDs.add(ID)
                    else:
                        #if it's not compatible with the new SNP, then it loses its second chance
                        no_more_chances_IDs.add(ID)
                
                #reinstate the second chancers
                for ID in second_chance_IDs:
                    clusters_forward[i]["IDs"].add(ID)
                    knockouts["IDs"].remove(ID)
                    clusters_forward[i]["nearest"][ID] = knockouts["nearest"].pop(ID)
                    clusters_forward[i]["n"][ID] = knockouts["n"].pop(ID)
                    del knockouts["by"][ID]
                
                #delete the no more chancers
                for ID in no_more_chances_IDs:
                    knockouts["IDs"].remove(ID)
                    del knockouts["nearest"][ID]
                    del knockouts["n"][ID]
                    del knockouts["by"][ID]
                
                #finally, record the new knockout IDs
                for ID in new_knockout_IDs:
                    knockouts["IDs"].add(ID)
                    knockouts["nearest"][ID] = clusters_forward[i-1]["nearest"][ID]
                    knockouts["n"][ID] = clusters_forward[i-1]["n"][ID]
                    knockouts["by"][ID] = matches[i]
    
    if second_chances: knockouts = {"IDs":set(), "nearest":{}, "n":{}, "by":{}}

    if not silent: print("\n    Backward search for compatible SNPs\n    ", file=sys.stderr)
    
    for i in reversed(range(0, total-1)):
        if not silent and i % onePercent == 0: print(".", end="", file=sys.stderr, flush=True)
        
        #if this SNP is in the previous cluster, import everything
        if matches[i] in clusters_reverse[i+1]["IDs"]:
            clusters_reverse[i]["IDs"].update(clusters_reverse[i+1]["IDs"])
            clusters_reverse[i]["nearest"].update(clusters_reverse[i+1]["nearest"])
            clusters_reverse[i]["nearest"][matches[i]] = positions[i]
            clusters_reverse[i]["n"].update(clusters_reverse[i+1]["n"])
            clusters_reverse[i]["n"][matches[i]] += 1
        
        #if current SNP pattern is not in the previous cluster, only import compatibles and send knockouts to the knockout dictionary
        else:
            if second_chances: new_knockout_IDs = set()
            for ID in clusters_reverse[i+1]["IDs"]:
                if passes_derived_gamete_test(patterns[ID], patterns[matches[i]], ploidies):
                    clusters_reverse[i]["IDs"].add(ID)
                    clusters_reverse[i]["nearest"][ID] = clusters_reverse[i+1]["nearest"][ID]
                    clusters_reverse[i]["n"][ID] = clusters_reverse[i+1]["n"][ID]
                elif second_chances:
                    new_knockout_IDs.add(ID)
            
            if second_chances:
                #check the knockouts in case some can be brought back
                second_chance_IDs = set()
                no_more_chances_IDs = set()
                for ID in knockouts["IDs"]:
                    #first see if it's compatible with the new SNP
                    if passes_derived_gamete_test(patterns[ID], patterns[matches[i]], ploidies):
                        #now see if its nemesis has been knocked out
                        if knockouts["by"][ID] in new_knockout_IDs:
                            #if we get here, this knockout should be compatible with the new cluster, so give it a second chance
                            second_chance_IDs.add(ID)
                    else:
                        #if it's not compatible with the new SNP, then it loses its second chance
                        no_more_chances_IDs.add(ID)
                
                #reinstate the second chancers
                for ID in second_chance_IDs:
                    clusters_reverse[i]["IDs"].add(ID)
                    knockouts["IDs"].remove(ID)
                    clusters_reverse[i]["nearest"][ID] = knockouts["nearest"].pop(ID)
                    clusters_reverse[i]["n"][ID] = knockouts["n"].pop(ID)
                    del knockouts["by"][ID]
                
                #delete the no more chancers
                for ID in no_more_chances_IDs:
                    knockouts["IDs"].remove(ID)
                    del knockouts["nearest"][ID]
                    del knockouts["n"][ID]
                    del knockouts["by"][ID]
                
                #finally, record the new knockout IDs
                for ID in new_knockout_IDs:
                    knockouts["IDs"].add(ID)
                    knockouts["nearest"][ID] = clusters_reverse[i+1]["nearest"][ID]
                    knockouts["n"][ID] = clusters_reverse[i+1]["n"][ID]
                    knockouts["by"][ID] = matches[i]
    
    if not silent: print("\n    Reconciling left and right clusters\n    ", file=sys.stderr)
    
    clusters_combined = [{"IDs":set(), "n":{}} for i in range(total)]
    
    #reconcile forward and reverse by iterating from nearest first
    for i in range(total):
        
        if not silent and i % onePercent == 0: print(".", end="", file=sys.stderr, flush=True)
        
        forward_unique_IDs = clusters_forward[i]["IDs"] - clusters_reverse[i]["IDs"]
        reverse_unique_IDs = clusters_reverse[i]["IDs"] - clusters_forward[i]["IDs"]
        #if one of these is empty, just combine all (because one is a subset of the other)
        if len(reverse_unique_IDs) == 0:
            clusters_combined[i]["IDs"].update(clusters_forward[i]["IDs"])
            clusters_combined[i]["n"].update(clusters_forward[i]["n"])
            for ID in clusters_reverse[i]["IDs"]:
                clusters_combined[i]["n"][ID] += clusters_forward[i]["n"][ID]
        elif len(forward_unique_IDs) == 0:
            clusters_combined[i]["IDs"].update(clusters_reverse[i]["IDs"])
            clusters_combined[i]["n"].update(clusters_reverse[i]["n"])
            for ID in clusters_forward[i]["IDs"]:
                clusters_combined[i]["n"][ID] += clusters_reverse[i]["n"][ID]
        else:
            #otherwise get intersection and add closest patternID iteratively
            clusters_combined[i]["IDs"].update(clusters_forward[i]["IDs"].intersection(clusters_reverse[i]["IDs"]))
            for ID in clusters_combined[i]["IDs"]:
                clusters_combined[i]["n"][ID] = clusters_forward[i]["n"][ID] + clusters_reverse[i]["n"][ID]
            
            #Now, for the remaining mutations, we get the distance to the nearest mutation of each type
            left_IDs = np.array(list(forward_unique_IDs))
            left_distances = np.array([positions[i] - clusters_forward[i]["nearest"][ID] for ID in left_IDs])
            left_order = np.argsort(left_distances)
            left_IDs = left_IDs[left_order]
            left_distances = left_distances[left_order]
            
            right_IDs = np.array(list(reverse_unique_IDs))
            right_distances = np.array([clusters_reverse[i]["nearest"][ID] - positions[i] for ID in right_IDs])
            right_order = np.argsort(right_distances)
            right_IDs = right_IDs[right_order]
            right_distances = right_distances[right_order]
            
            #now add them iteratively
            left_IDs_to_add = []
            right_IDs_to_add = []
            i_left = 0
            max_left = len(left_IDs)
            i_right = 0
            max_right = len(right_IDs)
            while True:
                if i_left == max_left:
                    if i_right == max_right:
                        break
                    else:
                        if compatible_with_all(right_IDs[i_right], left_IDs_to_add, patterns, ploidies):
                            right_IDs_to_add.append(right_IDs[i_right])
                        i_right += 1
                elif i_right == max_right:
                    if compatible_with_all(left_IDs[i_left], right_IDs_to_add, patterns, ploidies):
                        left_IDs_to_add.append(left_IDs[i_left])
                    i_left += 1
                else:
                    #if we get here, neither are maxed, so check which is closer
                    if left_distances[i_left] <= right_distances[i_right]:
                        if compatible_with_all(left_IDs[i_left], right_IDs_to_add, patterns, ploidies):
                            left_IDs_to_add.append(left_IDs[i_left])
                        i_left += 1
                    else:
                        if compatible_with_all(right_IDs[i_right], left_IDs_to_add, patterns, ploidies):
                            right_IDs_to_add.append(right_IDs[i_right])
                        i_right += 1
            #finally, add these to the set
            clusters_combined[i]["IDs"].update(set(left_IDs_to_add + right_IDs_to_add))
            for ID in left_IDs_to_add: clusters_combined[i]["n"][ID] = clusters_forward[i]["n"][ID]
            for ID in right_IDs_to_add: clusters_combined[i]["n"][ID] = clusters_reverse[i]["n"][ID]
        
        clusters_combined[i]["n"][matches[i]] -= 1 #subtract 1 for the focal ID so as not to count it twice
   
    del clusters_forward
    del clusters_reverse
    
    if not silent: print("\n    Merging adjacent identical intervals.\n    ", file=sys.stderr)
    #to make tree intervals, merge adjacent intervals with the same patternIDs
    clusters_merged = [clusters_combined[0]]
    clusters_merged[0]["interval_idx"] = [0, 0]
    
    for i in range(1, total):
        
        if not silent and i % onePercent == 0: print(".", end="", file=sys.stderr, flush=True)
        
        if clusters_merged[-1]["IDs"] == clusters_combined[i]["IDs"]:
            clusters_merged[-1]["interval_idx"][1] = i
        else:
            clusters_merged.append(clusters_combined[i])
            clusters_merged[-1]["interval_idx"] = [i, i]
    
    n_trees = len(clusters_merged)
    
    for i in range(n_trees):
        clusters_merged[i]["interval"] = [positions[j] for j in clusters_merged[i]["interval_idx"]]
        clusters_merged[i]["interval_ext"] = clusters_merged[i]["interval"][:]
    
    if seq_start: clusters_merged[0]["interval_ext"][0] = seq_start
    if seq_len: clusters_merged[-1]["interval_ext"][1] = seq_len
    
    #extend intervals to remove gaps
    for i in range(1, n_trees):
        new_start = np.ceil((clusters_merged[i-1]["interval"][1] + clusters_merged[i]["interval"][0])/2) 
        clusters_merged[i]["interval_ext"][0] = new_start
        clusters_merged[i-1]["interval_ext"][1] = new_start-1
    
    if not silent: print(f"\n    {n_trees} tree intervals identified.\n    ", file=sys.stderr)
    
    return clusters_merged


class Tree:
    def __init__(self, n_leaves, interval=None):
        self.leaves = list(range(n_leaves)) #integers from zero
        self.root = n_leaves # root is always the very next integer
        self.parents = [self.root] #a list of all non-leaf nodes
        self.node_children = {self.root: self.leaves} #start with a polytomy
        for leaf in self.leaves: self.node_children[leaf] = []
        self.node_parent = dict([(leaf, self.root) for leaf in self.leaves])
        self.next_node = self.root + 1
        self.interval = interval
    
    def children(self, node): #just added this so that the tree class behaves like a tskit tree
        return self.node_children[node]
    
    def add_node(self, children):
        new_node = self.next_node
        #this will add a node above the listed children and below their mutual parent
        parent = self.node_parent[children[0]]
        #the children must already have the same parent. If hey do not, the node is no compaible with this tree
        for child in children[1:]: assert self.node_parent[child] == parent, "All children must have the same parent to begin with. Consider first using Tree.remove_node()."
        #add node with its new children
        self.node_children[new_node] = list(children)
        #remove these children from their parent node and add the new node instead
        self.node_children[parent] = [new_node] + [c for c in self.node_children[parent] if c not in children]
        #add parent for the new node
        self.node_parent[new_node] = parent
        #correct parent for each child
        for child in children: self.node_parent[child] = new_node
        #add to parents
        self.parents.insert(0, new_node)
        #advence next node
        self.next_node += 1
        return new_node #just return its number in case this is useful
    
    def remove_node(self, node):
        #this will remove a node and attach its children to the node's parent
        new_parent = self.node_parent.pop(node)
        children = self.node_children.pop(node)
        #remove this child and add new add children to new parent
        self.node_children[new_parent].pop(self.node_children[new_parent].index(node))
        self.node_children[new_parent] += children
        for child in children: self.node_parent[child] = new_parent
        self.parents.pop(self.parents.index(node))
    
    def _newick_(self, node=None, node_labels=None):
        if not node: node = self.root
        _string = "("
        i=0
        for child in self.node_children[node]:
            if i > 0: _string += "," 
            if self.node_children[child] == []:
                _string += str(node_labels[child] if node_labels and child in node_labels else child)
            else:
                _string += self._newick_(node=child, node_labels=node_labels)
            i += 1
        
        return _string + ")"
    
    def as_newick(self, node_labels=None):
        #node_labels is a dictionary connecting numric node names to strings of your choice
        return self._newick_(node_labels=node_labels) + ";"
    
    def get_parent_child_edges(self, node_labels=None):
        #node_labels is a dictionary connecting numric node names to values of your choice
        edges = []
        for parent in self.parents:
            for child in self.node_children[parent]:
                if node_labels:
                    _parent_ = node_labels[parent] if parent in node_labels else parent
                    _child_ = node_labels[child] if child in node_labels else child
                    edges.append((_parent_,_child_,))
                else:
                    edges.append((parent,child,))
        
        return edges



def pick_children(desired_pattern, available_children_IDs, children_patterns, ignore_perfect_children=True):
    #each child had a pattern - the number of derived alleles in each tip
    #So the pattern for each haploid tip will be all zeros except a 1 where that tip is
    #patterns for deeper nodes have 1s in multiple slots
    #This function adds children until the sum of children patterns equals the desired one for the node
    picked_children = []
    picked_children_pattern_sum = np.zeros(len(desired_pattern), dtype=int)
    #children are considered in the order in which they are provided
    #remove any impossible children
    available_children_IDs = [ID for ID in available_children_IDs if (children_patterns[ID] - desired_pattern).max() == 0]
    #Perfect children are those who match the desired pattern exactly
    if ignore_perfect_children:
        available_children_IDs = [ID for ID in available_children_IDs if (children_patterns[ID] - desired_pattern).min() < 0]
    
    for ID in available_children_IDs:
        #add child if adding t doesn't exceed the number of derived alleles
        if (picked_children_pattern_sum + children_patterns[ID] - desired_pattern).max() <= 0:
            picked_children.append(ID)
            picked_children_pattern_sum += children_patterns[ID]
            if (picked_children_pattern_sum - desired_pattern).min() == 0: return picked_children
    
    return None



def patterns_to_tree(patterns_array, pattern_sums=None, pattern_indices=None, ploidies=None,
                     tree = None, node_pattern_idx_dict=None, node_pattern_dict=None,
                     multi_pass=True, skip_sanity_checks=True, return_pattern_dicts=False):
    
    if ploidies is None: ploidies = [1]*patterns_array.shape[1]
    if pattern_indices is None: pattern_indices = list(range(patterns_array.shape[0]))
    if pattern_sums is None: pattern_sums = patterns_array.sum(axis=1)
    
    n_ind = len(ploidies)
    #each mutation pattern gives the number of derived alleles per individual
    if not skip_sanity_checks:
        assert np.all(np.array([len(pat) for pat in patterns_array]) == n_ind), "All patterns must have the same number of individuals as the ploidies list"
        assert np.unique(patterns_array, axis=0).shape == patterns_array.shape, "There should be no repeats in the patterns"
        assert patterns_array.sum(axis=1).min() > 1, "Patterns must have at least two derived alleles"
        assert np.all(patterns_array <= ploidies), "No. derived alleles per individual cannot exceed individual ploidy"
    #node patterns give the number of derived mutations you would see per individual assuming a mutation at that node
    #the indices point to these node patterns
    indices_of_patterns_to_add = pattern_indices[:]
    #If there is no tree to begin with, we sart an empty one, as well as empty dict of pattern indices for each node and patterns for each node
    if not tree:
        n_leaves = sum(ploidies)
        tree = Tree(n_leaves=n_leaves)
        node_pattern_idx_dict = {}
        node_pattern_dict = {tree.root: ploidies}
        
        #Create the haploid node patterns for the leaf nodes
        node = 0
        for i in range(n_ind):
            for j in range(ploidies[i]): #one leaf for each ploidy level
                pattern = np.zeros(n_ind, dtype=int) #all zeros exept with a 1 at the position of the leaf
                pattern[i] = 1
                node_pattern_dict[node] = pattern
                node += 1
    
    else:
        #First remove nodes that will not appear in present tree 
        parents = [n for n in tree.parents if n != tree.root]
        for node in parents:
            idx = node_pattern_idx_dict[node]
            if idx not in pattern_indices:
                #this node should not be in the present tree
                tree.remove_node(node)
                node_pattern_idx_dict.pop(node)
                node_pattern_dict.pop(node)
            else:
                #Node(s) with this index are already in the tree, so no need to add it again
                try: indices_of_patterns_to_add.pop(indices_of_patterns_to_add.index(idx))
                except: pass
    
    #now sort the indices of patterns to add based on the total number of leaves under them, so that we traverse the tree downwards
    sums_of_patterns_to_add = pattern_sums[indices_of_patterns_to_add]
    indices_of_patterns_to_add = [indices_of_patterns_to_add[i] for i in np.argsort(sums_of_patterns_to_add)]
    
    #make a node and pick children for each pattern
    #Keep looping until there are no more free children to join into appropriate nodes
    change_made=True
    while change_made:
        change_made = False
        for idx in indices_of_patterns_to_add:
            desired_pattern = patterns_array[idx]
            picked_children = None
            for node in tree.parents:
                #only consider parent if its pattern is inclusive of the desired pattern
                pattern_difference = node_pattern_dict[node] - desired_pattern
                if np.any(pattern_difference < 0) or pattern_difference.sum() == 0: continue
                #only consider nodes with more than two children
                if len(tree.node_children[node]) < 3: continue
                #if we get here, we can try to pick children
                picked_children = pick_children(desired_pattern=desired_pattern, available_children_IDs=tree.node_children[node], children_patterns=node_pattern_dict)
                #if none, or all children of the node were picked, then there is nothing to do
                if not picked_children: continue
                #if we get here, we have viable children, so break
                break
            
            if picked_children:
                #new parent node
                new_node = tree.add_node(picked_children)
                node_pattern_idx_dict[new_node] = idx
                node_pattern_dict[new_node] = desired_pattern
                
                change_made = True
            
            else:
                if multi_pass: continue
                else: raise Exception(f"Mutation pattern {desired_pattern} not compatible with current tree.")
        
        if not multi_pass: break
    
    if return_pattern_dicts: return (tree, node_pattern_idx_dict, node_pattern_dict,)
    return tree


def infer_trees(patterns, ploidies, clusters, multi_pass=True, leaf_names=None, include_pattern_indices=False, silent=False):
    node_labels = dict(enumerate(leaf_names)) if leaf_names else None
    pattern_sums = patterns.sum(axis=1)
    total = len(clusters)
    onePercent = int(np.ceil(total/100))
    tree = None
    node_pattern_idx_dict = None
    node_pattern_dict = None
    for i in range(total):
        if not silent and i % onePercent == 0: print(".", end="", file=sys.stderr, flush=True)
        
        pattern_indices = list(clusters[i]["IDs"])
        
        tree, node_pattern_idx_dict, node_pattern_dict = patterns_to_tree(patterns_array=patterns, pattern_sums=pattern_sums,
                                                                          pattern_indices=pattern_indices, ploidies=ploidies,
                                                                          tree = tree, node_pattern_idx_dict = node_pattern_idx_dict,
                                                                          node_pattern_dict = node_pattern_dict, multi_pass=multi_pass,
                                                                          return_pattern_dicts=True)
        
        tree.interval = clusters[i]["interval_ext"]
        
        if include_pattern_indices: yield (tree, pattern_indices,)
        else: yield tree


def infer_ts(patterns, ploidies, clusters, multi_pass=True, silent=False):
    pattern_sums = patterns.sum(axis=1)
    total = len(clusters)
    n_leaves = sum(ploidies)
    tree = None
    node_pattern_idx_dict = None
    node_pattern_dict = None
    tableCollection = tskit.TableCollection(clusters[-1]["interval_ext"][1])
    for i in range(n_leaves):
        tableCollection.nodes.add_row(flags=1, time=0) #add node for each leaf and flag these as samples
    last_node = tableCollection.nodes.add_row(time=n_leaves) #Add root. Time is defined as number of descenents. Crude but works for now.
    current_edges = set()
    current_edge_left_pos = {}
    for cluster in clusters:
        
        pattern_indices = list(cluster["IDs"])
                
        tree, node_pattern_idx_dict, node_pattern_dict = patterns_to_tree(patterns_array=patterns, pattern_sums=pattern_sums,
                                                                          pattern_indices=pattern_indices, ploidies=ploidies,
                                                                          tree = tree, node_pattern_idx_dict = node_pattern_idx_dict,
                                                                          node_pattern_dict = node_pattern_dict, multi_pass=multi_pass,
                                                                          return_pattern_dicts=True)
        
        #add any new nodes to the node table
        for node in sorted([i for i in node_pattern_idx_dict.keys() if i > last_node]):
            _node_ = tableCollection.nodes.add_row(time=node_pattern_dict[node].sum()) #time is number of descendents. Crude but it works.
            assert node == _node_, f"Expected next node to be {_node_} but found {node} instead."
            last_node = node
        
        #add any edges if they are no longer in this tree
        parent_child_edges = set(tree.get_parent_child_edges())
        for edge in current_edges - parent_child_edges:
            #these ones are finished, so add them to the table
            tableCollection.edges.add_row(left=current_edge_left_pos.pop(edge), right=last_tree_right_pos,
                                          parent=edge[0], child=edge[1])
        
        for edge in parent_child_edges - current_edges:
            #these are new starting edges, so add left pos here
            current_edge_left_pos[edge] = cluster["interval_ext"][0]-1 #treesequence tables use 0-based positions, left inclusive right exclusive
        
        #update variables before moving on to next tree
        current_edges = parent_child_edges
        last_tree_right_pos = cluster["interval_ext"][1]
    
    #add remaining edges from the final tree
    for edge in current_edges:
        tableCollection.edges.add_row(left=current_edge_left_pos.pop(edge), right=last_tree_right_pos,
                                      parent=edge[0], child=edge[1])
    
    tableCollection.sort()
    
    return tableCollection.tree_sequence()


def makeHaploidNames(name,ploidy):
    if ploidy == 1: return([name])
    return [name + "_" + letter for letter in string.ascii_uppercase[:ploidy]]


def parsePloidyArgs(args, sampleIDs, default=2):
    if args.ploidy is not None:
        ploidies = [args.ploidy]*len(sampleIDs)
        ploidyDict = dict(zip(sampleIDs, ploidies))
    elif args.ploidy_file is not None:
        with open(args.ploidy_file, "rt") as pf:
            ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
            ploidies = [ploidyDict[ID] for ID in sampleIDs]
    else:
        ploidies = [default]*len(sampleIDs)
        ploidyDict = dict(zip(sampleIDs, ploidies))
    
    return(ploidies, ploidyDict,)

def main():
    ### parse arguments

    parser = argparse.ArgumentParser(prog="sticcs")
    
    subparsers = parser.add_subparsers(title = "subcommands", dest="prep_or_ts")
    subparsers.required = True
    
    prep_parser = subparsers.add_parser("prep", help="Prepare a vcf file with derive allele counts")
    
    prep_parser.add_argument("-i", "--input_vcf", help="Input vcf file", action = "store", default="-")
    prep_parser.add_argument("-o", "--output_vcf", help="Output vcf file", action = "store", default="-")
    prep_parser.add_argument("--outgroups", help="sample name(s) for outgroup(s)", action='store', nargs="+")
    prep_parser.add_argument("--use_REF_as_anc", help="Use REF allele as ancestral.", action = "store_true")
    prep_parser.add_argument("--samples", help="SampleIDs to include (NOTE, order is ignored)", action = "store", nargs = "+")
    
    ts_parser = subparsers.add_parser("ts", help="Make tree sequence")  
    
    ts_parser.add_argument("-i", "--input_vcf", help="Input vcf file made by sticcs prep", action = "store", default="-")
    ts_parser.add_argument("-o", "--out_prefix", help="Output file prefix", action = "store", required=True)
    ts_parser.add_argument("--output_format", help="Type of tree sequence output file(s)", choices=["newick", "tskit"], action = "store", default = "newick")
    ts_parser.add_argument("--output_mutation_info", help="Include information about mutation patterns in output", action='store_true')
    ts_parser.add_argument("--samples", help="SampleIDs to include (NOTE, order is ignored)", action = "store", nargs = "+")
    ts_parser.add_argument("--ploidy", help="Sample ploidy if all the same. Use --ploidy_file if samples differ.", action = "store", type=int)
    ts_parser.add_argument("--ploidy_file", help="File with samples names and ploidy as columns", action = "store")
    ts_parser.add_argument("--allow_second_chances", help="Consider SNPs that are separated by incompatible SNPs", action='store_true')
    ts_parser.add_argument("--single_pass", help="Single pass when building trees (only relevant for ploidy > 1)", action='store_true')
    ts_parser.add_argument("--silent", help="Do not report on progress", action='store_true')
    
    args = parser.parse_args()
    
    print(not sys.stdin.isatty(), file=sys.stderr)
    
    vcf = cyvcf2.VCF(args.input_vcf, samples=args.samples)
    
    if args.samples:
        for sampleID in args.samples: assert sampleID in vcf.samples, f"Sample ID {sampleID} not in vcf file."
    
    if args.prep_or_ts=="prep":
        print("\nMaking VCF file with DC field...", file=sys.stderr)
        make_vcf_with_DC_field(vcf, args.output_vcf, outgroups=args.outgroups, use_REF_as_anc=args.use_REF_as_anc)
        print("\nDone.", file=sys.stderr)
        exit()
    
    if args.output_format == "tskit":
        import tskit
    
    sampleIDs = vcf.samples
    
    ploidies, ploidyDict = parsePloidyArgs(args, sampleIDs)
    
    total_haps = sum(ploidies)
    
    dac_generator = parse_vcf_with_DC_field(vcf)
    
    leaf_names = [hapID for ID in sampleIDs for hapID in makeHaploidNames(ID, ploidyDict[ID])]
    
    chromLenDict = dict(zip(vcf.seqnames, vcf.seqlens))
    
    print(f"\nReading first chromosome...", file=sys.stderr)
    
    for chrom, positions, der_counts in dac_generator:
        
        assert positions.max() <= chromLenDict[chrom], f"\tVariant at position {positions.max()} exceeds chromosome length {chromLenDict[chrom]} for {chrom}."
        
        no_missing = np.all(der_counts >= 0, axis=1)
        
        site_sum = der_counts.sum(axis=1)
        informative = (1 < site_sum) & (site_sum < total_haps)
        
        usable_sites = np.where(no_missing & informative)
        
        #der counts and positions
        der_counts = der_counts[usable_sites]
        
        positions = positions[usable_sites]
        
        print(f"\nAnalysing {chrom}. {positions.shape[0]} usable SNPs found.", file=sys.stderr)
        
        patterns, matches, n_matches  = get_patterns_and_matches(der_counts)
        
        print("\nFinding clusters...", file=sys.stderr)
        clusters = get_clusters(patterns, matches, positions, ploidies=ploidies,
                                second_chances=args.allow_second_chances,
                                seq_start=1, seq_len=chromLenDict[chrom],
                                silent=args.silent)
        
        print("\nBuilding trees and writing output...", file=sys.stderr)
        
        if args.output_format == "tskit":
            ts = infer_ts(patterns, ploidies, clusters, multi_pass = not args.single_pass, silent=args.silent)
            ts.dump(args.out_prefix + "." + chrom + ".ts")
        
        else:
            tree_generator = infer_trees(patterns, ploidies, clusters, multi_pass = not args.single_pass, leaf_names=leaf_names,
                                        include_pattern_indices=True, silent=args.silent)
            
            with open(args.out_prefix + "." + chrom + ".windowdata.tsv", "wt") as out_data:
                out_data.write("chrom\tstart\tend")
                if args.output_mutation_info: out_data.write("\tpattern_indices")
                out_data.write("\n")
                with gzip.open(args.out_prefix + "." + chrom + ".trees.gz", "wt") as out_trees:
                        for tree, pattern_indices in tree_generator:
                            out_trees.write(tree.as_newick(node_labels=leaf_names) + "\n")
                            out_data.write("\t".join([chrom, str(tree.interval[0]), str(tree.interval[1])]))
                            if args.output_mutation_info: out_data.write("\t" + ",".join([str(i) for i in pattern_indices]))
                            out_data.write("\n")
            
            if args.output_mutation_info:
                np.savetxt(args.out_prefix + "." + chrom + ".patterns.tsv.gz", np.column_stack([n_matches, patterns]), delimiter="\t", fmt='%i')
        
        print("\nDone\n", file=sys.stderr)
    


if __name__ == '__main__':
    main()
