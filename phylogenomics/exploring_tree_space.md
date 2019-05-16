Exploring the tree space to find a (near) optimal solution

![tree_space]

Tree space = distribution of trees in a 3 space - optimal solutions represented as peaks
Several islands of optimal/suboptimal solutions

Related trees - close together in the tree space (degree of relatedness)

It is not possible to move downwards -> so to change from peak to peak you need a different starting point
Always moves towards a better solution

# Mechanics of exploring the tree space:

## 1). Select a starting tree

Starting tree can be:
- a tree buid from distance : NJ tree
- a random tree

## 2). Searching the tree space

### Branches swapping techniques: Iterative techniques

From low to high efficiency and computing amount required
- `NNI: nearest neighboor interchange`:
  > each branch is considered as a subtree
  > 1 pair of branches swapped at each iteration
- SPR: random subtree prunning and regraphting
  > untill all the possible regraphting are done
- TBR: tree bissection and reconnection
  > cut on interior branch - reconnect at the end subtree
  > to a local maximum

Stops when no better (according to a [criteria] is optained)
Q ? for first one or all ?

[criteria]: Link to that

## 3). From local to "true" optimal peak

- randomizing starting tree: `Rachet` (Kevin Nixon)
  > Both for ML and Parsimony
  > pushes the search into different area of the tree space
  - randomly changes the weights ( = transformation parameters
`character evolution/changes matrix`)
  - 1 starting tree -> 2 randomly weighted charact changed trees
