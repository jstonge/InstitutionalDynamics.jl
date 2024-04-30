## InstitutionalDynamics.jl

We are using group-based master equation models (GMEs) to study how array of institutional responses arise from group dynamics.

## Installation

Clone this repo. Make sure you have `julia=1.10` and higher installed. If you need to update your `Julia`, see the [juliaup project](https://github.com/JuliaLang/juliaup). We prefer most recent version because precompile time are shorter.

## Usage

 ```shell
 # the rest of the params will be the default below
 julia models/source-sink1.jl --beta 0.27 -g 1.1 
 ```
 Running `julia models/source-sink1.jl --help` will give you the argument names and current default values:
 
 ```shell
 usage: sourcesink1.jl [--db DB] [-L L] [-O O] [--beta BETA] [-g G]
                      [-r R] [-b B] [-c C] [-m M] [-o O] [-h]

optional arguments:
  --db DB      Use Database to query parameters
  -L L         LIMIT of rows (type: Int64, default: 5)
  -O O         The OFFSET clause after LIMIT specifies how many rows
               to skip at the beginning of the result set. (type:
               Int64, default: 0)
  --beta BETA  Spreading rate from non-adopter to adopter beta (type:
               Float64, default: 0.07)
  -g G         Recovery rate gamma, i.e. rate at which adopters loose
               behavioral trait (type: Float64, default: 1.0)
  -r R         Global behavioral diffusion rho (allows the behaviour
               to spread between groups) (type: Float64, default: 0.1)
  -b B         Group benefits b (type: Float64, default: 0.18)
  -c C         Institutional cost c (type: Float64, default: 1.05)
  -m M         Noise u (type: Float64, default: 0.0001)
  -o O         Output file for results (default: ".")
  -h, --help   show this help message and exit
 ```
 Once this is done running, the rendered file should reflect the changes. To see the changes on the web page, we need to push the changes on Github. Then Github action will take care of update the web page.
