# PEF

## installation

* PEF has dependencies on [SDSL-library](https://github.com/Hiipivahalko/sdsl-lite.git) and [SUX](https://github.com/vigna/sux.git).
Install libraries to your include path (see instructions from links).
* All PEF implementations are template, to use PEF, copy the content from `./src/include` directory to your include path.


### Dev config

For VIM run following command after external lib installation for coc.vim to find sources:

```bash
bear -- make
```
This command creates file `compile_commands.json`.

## TEST libraries

### test setup

```bash
mkdir -p external/cpp_libs
cd external
git clone https://github.com/Hiipivahalko/sdsl-lite.git
cd sdsl-lite && ./install ../cpp_libs && cd .. # install sdsl to cpp_libs directory
git clone https://github.com/vigna/sux.git
cp -a sux/sux cpp_libs/include
```

### running tests

Run scriptfile `run_all_tests.sh` in project root directory:
```bash
./run_all_tests.sh
```

For running a specific testsets:
```bash
./run_all_tests.sh input <prefix_of_test_case>
```

## PEF Uniform

### rank operation

#### Description

Given an position `i`, do `rank(i)`

#### Algorithm

The rank operation first of all has to do a rank in the vector `L` to know in which block the bit i is.

Then, if the block is not the first block, it has to sum `b` ones times the number of blocks before the actual block, this is because before this block there is `b` ones for each block.

Then, depending of the value of `B[blk]` and `P[blk]` it will use different rank operations:
* if `B[blk] = 1` the block uses Elias-Fano Encoding, so we will use rank operation of that encoding, then:
  * if the block is the first block it can do the operation rank with the value of `i`
  * if the block is not the first block, then it has to know the position of the `i` bit in that block, to do that is has to do the difference between `i - 1`and the result of a select in `L` with the value of the block, to know what was the position of the last `1` from the previous block, and finally we can do rank with that value
* if `B[blk] = 0 && P[blk] != NULL` the block uses a plan bit vector encoding, so we will use rank operation of that encoding, then:
  * if the block is the first block it can do the operation rank with the value of `i`
  * if the block is not the first block, then is has to know the position of the `i` bit in that block, to do that it has to do the difference between `i - 1` and the result of a select in `L` with the value of the block, to know what was the position of the last `1` from the previous block, and finally we can do rank with that value
* if `B[blk] = 0 && P[blk] == NULL` the block is full of ones, then:
  * if the block is the first block the result of the rank is just `i`
  * if the block is not the first block, it has to know the position of the `i` bit in that block, to do that have to do the difference between `i - 1` and the result of a select in `L` with the value of the block, to know what was the position of the last `1` from the previous block, and the rank of in that block will be that value
