#!/bin/bash
# Status check for all BatteryMat prospective delithiation chains.
# Run from /data/jlee859/dft/neurips:  bash check_chains.sh
# Read-only: inspects dft_inputs/, energies.json, logs, and squeue. Changes nothing.

cd "$(dirname "$0")"
printf "%-14s %-14s %6s %8s %9s %-9s %s\n" \
    "JID" "latest_step" "nLi0" "recorded" "voltage" "queue" "status"
echo "--------------------------------------------------------------------------------"

DONE=0; ACTIVE=0; STALLED=0
for j in $(ls dft_inputs); do
    sup=$(ls -d dft_inputs/"$j"/supercell_* 2>/dev/null | head -1)
    [ -z "$sup" ] && continue

    latest=$(ls -d "$sup"/step_* 2>/dev/null | sort -t_ -k2 -n | tail -1 | xargs basename)
    nli0=$(ls -d "$sup"/step_00_Li* 2>/dev/null | head -1 | sed 's/.*_Li//')
    nrec=$(python3 -c "import json;print(len(json.load(open('$sup/energies.json'))))" 2>/dev/null || echo 0)
    [ -f "$sup/voltage_curve.png" ] && volt="yes" || volt="no"

    q=$(squeue -h -u "$USER" -n "bm_$j" -o "%T" 2>/dev/null | head -1)
    [ -z "$q" ] && q="-"

    # last chain log for this jid, check its error stream
    lastlog=$(ls -t chain_bm_"$j"_*.err 2>/dev/null | head -1)
    errsig=""
    [ -n "$lastlog" ] && errsig=$(grep -l -iE "segfault|error|abort|traceback|CANCELLED|time limit" "$lastlog" 2>/dev/null)

    if [ "$volt" = "yes" ]; then
        st="COMPLETE"; DONE=$((DONE+1))
    elif [ "$q" != "-" ]; then
        st="active ($q)"; ACTIVE=$((ACTIVE+1))
    elif [ -n "$errsig" ]; then
        st="STALLED - last log has errors: $lastlog"; STALLED=$((STALLED+1))
    else
        st="STALLED - no job queued, no voltage curve"; STALLED=$((STALLED+1))
    fi

    printf "%-14s %-14s %6s %8s %9s %-9s %s\n" "$j" "$latest" "$nli0" "$nrec" "$volt" "$q" "$st"
done

echo "--------------------------------------------------------------------------------"
echo "COMPLETE: $DONE   active: $ACTIVE   STALLED: $STALLED   (of $(ls dft_inputs | wc -l))"
echo
echo "Final voltages from completed chains (from most recent logs):"
for j in $(ls dft_inputs); do
    log=$(grep -l "V=" chain_bm_"$j"_*.out 2>/dev/null | tail -1)
    [ -n "$log" ] && echo "  $j: $(grep -E "x=.*V=" "$log" | tail -3 | tr '\n' ' ')"
done
echo
if [ "$STALLED" -gt 0 ]; then
    echo "For each STALLED jid: check its last .err, then resubmit with"
    echo "  sbatch --job-name=bm_<JID> chain_step.sh <JID>"
    echo "(completed steps are skipped automatically)."
fi