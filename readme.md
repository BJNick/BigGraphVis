GraphVis test implementation uses Brinkmann implementation (cited below).

* G.G. Brinkmann, [K.F.D. Rietveld](https://liacs.leidenuniv.nl/~rietveldkfd) and [F.W. Takes](https://liacs.leidenuniv.nl/~takesfw), [Exploiting GPUs for fast force-directed visualization of large-scale networks](https://dx.doi.org/10.1109/ICPP.2017.47), in Proceedings of the 46th International Conference on Parallel Processing (ICPP), pp. 382-391, 2017.


#### System Requirements

A CUDA capable GPU. Currently only Linux is supported.
The windows Desktop application under developing and it will be published here.

#### Obtaining all code
This repository contains a submodule (`lib/pngwriter`). Be sure to run
```
git submodule init && git submodule update
```
from the root of this Git repository before compiling. The code also depends on the `libpng` library (including its development headers). It should be possible to obtain this using the package manager for your Linux distribution.

#### Compiling
A `Makefile` is located in `builds/linux`. Running
```
v
```
from this directory compiles `graph_viewer` with CUDA support.
To compile without CUDA support, run `make graph_viewer CUDA_SUPPORT=0`.

#### Usage
`graph_viewer gpu|cpu max_iterations num_snaps sg|wg scale gravity exact|approximate edgelist_path out_path [png|csv|bin]`

`gpu`            : choose a parallel GPU implementation .

`max_iterations`     : how many iterations of the layout algorithm to run

`num_snaps`          : choose how many times during the layout process a visualization should be rendered

`wg|sg`              : choose between weak gravity (inversely proportional to distance) or
                     strong gravity

`scale`              : scale repulsive force

`gravity`            : scale gravitational force

`exact|approximate`  : choose between the exact/pairwise O(|V|^2) repulsive force calculation or the O(|V|lg(|V|))
                     approximation using Barnes-Hut (GPU implementation only supports Barnes-Hut)

`edgelist_path`      : Text file (ascii) containing node IDs for each edge on a separate line (whitespace separated).
                       Lines starting with a `#`, the direction of edges, and self-loops are ignored.

`out_path`           : path to write resulting layout to

`[png|csv|bin]` is optional, defaulting to `png`, and determines the format of the layout written to `out_path`.

#### References
<a name="jacomy14"><sup>1</sup></a> M. Jacomy, T. Venturini, S. Heymann, and M. Bastian, ["Forceatlas2, a continuous graph layout algorithm for handy network visualization designed for the Gephi software"](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679), PLoS ONE, vol. 9, no. 6, pp. 1–12, 2014.

<a name="bastian09"><sup>2</sup></a> M. Bastian, S. Heymann, and M. Jacomy, ["Gephi: an open source software for exploring and manipulating networks."](https://aaai.org/ocs/index.php/ICWSM/09/paper/view/154) in Proceedings of International Conference on Web and Social Media (ICWSM), 2009, pp. 361–362.

<a name="barnes86"><sup>3</sup></a>J. Barnes and P. Hut, ["A hierarchical O(N log N) force-calculation algorithm"](https://www.nature.com/nature/journal/v324/n6096/abs/324446a0.html), Nature, vol. 324, pp. 446–449, 1986.

<a name="burtscher11"><sup>4</sup></a> M. Burtscher and K. Pingali, ["An efficient CUDA implementation of the tree-based Barnes Hut n-body algorithm"](https://www.sciencedirect.com/science/article/pii/B9780123849885000061), in GPU Computing Gems Emerald Edition, W. mei W. Hwu, Ed., 2011, ch. 6, pp. 75–92.

#### License
Most source files for this program are released under the GNU Affero General Public License. The license notice in each file provides more information. A copy of the GNU Affero General Public License can be found in the `LICENCE` file.

#### Disclaimer
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
