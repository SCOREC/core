#!/bin/bash -ex
scan-build -o /lore/dibanez/cdash/static_analysis_results /usr/local/bin/gmake "$@"
