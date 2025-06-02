import sys
import sphinx.cmd.build

sys.exit(sphinx.cmd.build.main(["-b", "html", "docs/source", "docs/_build/html"]))
