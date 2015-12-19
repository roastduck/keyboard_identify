Keyboard Identify (OLD IO)
========

Given some iOS keyboard screenshots (PNG), output the positions of the keys.

Output : YAML Mode

Dependencies
--------
- boost-filesystem
- cpp-yaml

Compile
--------
run `make yamlio` to compile into keyboard

Run
--------
run `./keyboard` to parse folders in current directory.

run `./keyboard <path>` to parse folders in _path_, for example, `./keyboard examples`.
