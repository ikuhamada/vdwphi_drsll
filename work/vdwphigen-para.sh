#!/bin/sh

BINDIR=.
BIN=${BINDIR}/vdwphigen_drsll

MPIRUN=mpirun
NPROCS=08

LOG=log
ERR_LOG=err

if [ ! -e ${BIN} ];
then
  ${BIN} does not exist.
  exit
fi

# # for plotting
# dmax=30.d0
# mD=300
# delmax=0.5d0
# mdel=5
# ${MPIRUN} -np ${NPROCS} ${BIN} -Dmax ${dmax} -mD ${mD} -delmax ${delmax} -mdel ${mdel} -plot

# for testing
LOG=log_plot
dmax=10.d0
mD=100
delmax=0.9d0
mdel=9
(${MPIRUN} -np ${NPROCS} ${BIN} -Dmax ${dmax} -mD ${mD} -delmax ${delmax} -mdel ${mdel} > ${LOG} -plot ) >& ${ERR_LOG}

# for production
LOG=log
( ${MPIRUN} -np ${NPROCS} ${BIN} > ${LOG} ) >& ${ERR_LOG}

