# MIT 6.034 Lab 8: Bayesian Inference
# Written by Dylan Holmes (dxh), Jessica Noss (jmn), and 6.034 staff

from nets import *


#### ANCESTORS, DESCENDANTS, AND NON-DESCENDANTS ###############################

def get_ancestors(net, var):
    "Return a set containing the ancestors of var"
    ancestors = set()
    for ancestor in net.get_parents(var):
        ancestors.add(ancestor)
        ancestors = ancestors.union(get_ancestors(net,ancestor))
    return ancestors

def get_descendants(net, var):
    "Returns a set containing the descendants of var"
    descendants = set()
    for descendant in net.get_children(var):
        descendants.add(descendant)
        descendants = descendants.union(get_descendants(net,descendant))
    return descendants

def get_nondescendants(net, var):
    "Returns a set containing the non-descendants of var"
    descendants = get_descendants(net,var)
    descendants.add(var)
    return set(net.get_variables()) - descendants
    
    

def simplify_givens(net, var, givens):
    """If givens include every parent of var and no descendants, returns a
    simplified list of givens, keeping only parents.  Does not modify original
    givens.  Otherwise, if not all parents are given, or if a descendant is
    given, returns original givens."""
    if givens is None: # nothing to simplify here
        return givens
    simplified_givens = {}
    descendants = get_descendants(net,var)
    parents = net.get_parents(var)
    if bool(set(givens) & descendants) or not (parents.issubset(set(givens))):
        return givens
    bad_parents = set(givens) - parents
    for k,v in givens.items():
        if not k in bad_parents:
            simplified_givens[k] = v
    return simplified_givens


#### PROBABILITY ###############################################################

def probability_lookup(net, hypothesis, givens=None):
    "Looks up a probability in the Bayes net, or raises LookupError"
    givens = simplify_givens(net, list(hypothesis)[0], givens)
    try:
        return net.get_probability(hypothesis, givens)
    except:
        raise LookupError
    

def probability_joint(net, hypothesis):
    "Uses the chain rule to compute a joint probability"
    
    givens = {}
    probability = 1
    for var in net.topological_sort():
        if var in hypothesis:
            probability *= probability_lookup(net, {var: hypothesis[var]}, givens)
            givens[var] = hypothesis[var]
    return probability


def probability_marginal(net, hypothesis):
    "Computes a marginal probability as a sum of joint probabilities"
    probability = 0
    partitions = net.combinations(net.get_variables(), hypothesis)
    for partition in partitions:
        probability += probability_joint(net,partition)
    return probability

def probability_conditional(net, hypothesis, givens=None):
    "Computes a conditional probability as a ratio of marginal probabilities"  
    if givens is None:
        return probability_marginal(net,hypothesis)
    
    for var in hypothesis:
        if var in givens:
            if hypothesis[var] != givens[var]:
                return 0

    denominator = probability_marginal(net,givens)         
    givens = dict(givens, **hypothesis)
    numerator = probability_marginal(net,givens)
    
    return numerator/denominator
    
def probability(net, hypothesis, givens=None):
    "Calls previous functions to compute any probability"
    return probability_conditional(net,hypothesis,givens)


#### PARAMETER-COUNTING AND INDEPENDENCE #######################################

def number_of_parameters(net):
    "Computes minimum number of parameters required for net"
    vars = net.get_variables()
    # represent each var as (domain of var, product[domain of each parents of var])
    vars_details = map(lambda var: (len(net.get_domain(var)),\
                                    product(map(lambda p: len(net.get_domain(p)),\
                                                net.get_parents(var)))),vars)
    vars_totals = map(lambda (dom,prod): (dom-1)*max(prod,1),vars_details)
    return sum(vars_totals)

                                    


def is_independent(net, var1, var2, givens=None):
    """Return True if var1, var2 are conditionally independent given givens,
    otherwise False.  Uses numerical independence."""
    # v1 & v2 independant given "givens" if P(v1,v2|givens) = P(v1|givens)*P(v2|givens)
    for domain1 in net.get_domain(var1):
        for domain2 in net.get_domain(var2):
            p_left = probability(net, {var1: domain1, var2: domain2}, givens)
            p_right = probability(net,{var1:domain1},givens)*\
                      probability(net,{var2:domain2},givens) 
            if not approx_equal(p_left,p_right, 0.0000000001):
                return False
    return True
                                  


    

def is_structurally_independent(net, var1, var2, givens=None):
    """Return True if var1, var2 are conditionally independent given givens,
    based on the structure of the Bayes net, otherwise False.
    Uses structural independence only (not numerical independence)."""
    # ancestralize
    vars = set([var1,var2])
    if givens is not None:
        vars.update(set(givens))
    ancestors = set()
    for var in vars:
        ancestors.update(get_ancestors(net,var))
    vars.update(ancestors)
    subnet = net.subnet(vars)
    
    # moralize
    for var in subnet.topological_sort():
        parents = list(subnet.get_parents(var))
        for i in range(len(parents)):
            for j in range(i, len(parents)):
                if i != j:
                    subnet.link(parents[i], parents[j])

    # disorient
    subnet.make_bidirectional()

    # delete givens & their edges
    if givens is not None:
        for var in givens:
            subnet.remove_variable(var)

    # check var1 and var2
    if subnet.find_path(var1,var2) is None:
        return True
    return False
        
    
#### SURVEY ####################################################################

NAME = None
COLLABORATORS = None
HOW_MANY_HOURS_THIS_LAB_TOOK = None
WHAT_I_FOUND_INTERESTING = None
WHAT_I_FOUND_BORING = None
SUGGESTIONS = None
