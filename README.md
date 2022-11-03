# Restriction and transcription motifs

### Task 1 - Double Digest Problem (DDP)
* Implement function for brute force DDP algorithm for one fragments arrangement.
* Next, rework the function to check all possible fragments arrangements.


### Task 2 - Partial Digest Problem (PDP)
Implement recursive function for PDP algorithm according to following pseudocode:

```
PartialDigestProblem(L)
1  width <- Maximum element in L
2  Delete(width,L)
3  X <- {0, width}
4  Place(L, X)
```

``` 
Place(L, X)
1   if L is empty
2     output X
3     return
4   y <- Maximum element in L 
5   if delta(y, X) is a subset of L
6     Add y to X and remove lengths delta(y, X) from L
7     Place(L, X)
8     Remove y from X and add lengths delta(y, X) to L
9   if delta(width - y, X) is a subset of L
10    Add width − y to X and remove lengths delta(width − y, X) from  L
11    Place(L, X)
12    Remove width − y from X and add lengths delta(width − y, X) to L
13  return
```
Note: Algorithm works with the list of pairwise distances, `L`, and uses the function `Delete(y, L)` which removes the value `y` from `L`. 
We use the notation `delta(y, X)` to denote the multiset of distances between a point `y` and all points in a set `X`.


### Task 3 - Brute force motif search
Implement brute force motif search algorithm.


1. Function `Score()`

   Inputs:
   * a DNAStringSet of sequences (for example file `seq_score.fasta`)
   * a list of initial indexes
   * a length of a motif
   
   Outputs:
   * the best score
   * a block of sequences

2. Function `NextLeaf()`
    ```NextLeaf(a, L, k)
    1   for i ← L to 1
    2       if ai < k
    3           ai ← ai + 1
    4           return a
    5       ai ← 1
    6   return a
    ```
    Inputs:
   * a - `rep(1,L)`
   * L - number of DNA sequences
   * k - `n - l + 1`
   * n - length of each DNA sequence
   * l - motif length

3. Function `BFMotifSearch()`
    ```
    1   s ← (1, 1, . . . , 1)
    2   bestScore ← Score(s, DNA)
    3   while forever
    4       s ← NextLeaf(s, t, n −l + 1)
    5       if Score(s, DNA) > bestScore
    6           bestScore ← Score(s, DNA)
    7           bestMotif ← (s1, s2, . . . , st)
    8       if s = (1, 1, . . . , 1)
    9           return bestMotif
    ```
    Inputs:
   * DNA – a DNAStringSet of DNA sequences (for example file `seq_motif.fasta`)
   * t – number of DNA sequences
   * n – length of each DNA sequence
   * l – motif length

4. Function `NextVertex()`
    ```
    NextVertex(a, i, L, k)
    1   if i < L
    2       ai+1 ← 1
    3       return (a, i + 1)
    4   else
    5       for j ← L to 1
    6           if aj < k
    7               aj ← aj + 1
    8               return (a, j)
    9   return (a, 0)
    ```

5. Function `SimpleMotifSearch()`
    ```
    SimpleMotifSearch(DNA, t, n, l)
    1   s ← (1, . . . , 1)
    2   bestScore ← 0
    3   i ← 1
    4   while i > 0
    5       if i < t
    6           (s, i) ← NextVertex(s, i, t, n − l + 1)
    7       else
    8           if Score(s, DNA) > bestScore
    9               bestScore ← Score(s, DNA)
    10              bestMotif ← (s1, s2, . . . , st)
    11          (s, i) ← NextVertex(s, i, t, n − l + 1)
    12  return bestMotif
    ```


<details>
<summary>Download files from GitHub</summary>
<details>
<summary>Git settings</summary>

> * Configure the Git editor
> ```bash
> git config --global core.editor notepad
> ```
> * Configure your name and email address
> ```bash
> git config --global user.name "Zuzana Nova"
> git config --global user.email z.nova@vut.cz
> ```
> * Check current settings
> ```bash
> git config --global --list
> ```
>
</details>

* Create a fork on your GitHub account. 
  On the GitHub page of this repository find a <kbd>Fork</kbd> button in the upper right corner.
  
* Cloned forked repository from your GitHub page to a folder in your computer:
```bash
git clone <fork repository address>
```
* In a local repository, set new remote for project repository:
```bash
git remote add upstream https://github.com/mpa-prg/exercise_05.git
```

### Send files to GitHub
Create a new commit and send new changes to your remote repository.
* Add file to a new commit.
```bash
git add <file_name>
```
* Create commit, enter commit message, save the file and close it.
```bash
git commit
```
* Send new commit to your GitHub repository.
```bash
git push origin master
```
</details>
