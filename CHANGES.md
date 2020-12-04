#### v0.14.2

* Changed SFT cutoff short flag to `-c` from `-S` since `-S` conflicts with `--no-scale`.
* Bugfixes with displaying default values in help pages.
* `seidr top` now handles ties properly
* Fixed a bug where `seidr top` would not print top N edges if N is large and
it encounters a high scoring edge early

#### v0.14.1

* Fixed critical bug where seidr tools would only print version
* Added version print when logger is initialized
