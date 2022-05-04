#!/bin/sh

BINDIR=.
BIN=${BINDIR}/vdwphigen_drsll

dmax=30.d0
mD=300

delmax=0.5d0
mdel=5

${BIN} -Dmax ${dmax} -mD ${mD} -delmax ${delmax} -mdel ${mdel} -plot

