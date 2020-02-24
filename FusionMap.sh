#!/usr/bin/env bash
set -eo pipefail

usage() { echo "Usage: $0 [-t <threads>] [-o <outdir>] [-c <config>] [-g <genelist>] <reads> " 1>&2; exit 1; }

while getopts ":c:g:o:t:h" arg; do
    case "${arg}" in
        c)
            CONFIG="${OPTARG}"
            ;;
        g)
            GENELIST="${OPTARG}"
            ;;
        o)
            OUTPUTFOLDER="${OPTARG}"
            ;;
        t)
            THREADS="${OPTARG}"
            ;;
        h)
            usage
            exit 1
            ;;
        *)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 2
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${THREADS}" ] || [ -z "${OUTPUTFOLDER}" ] || [ -z "${CONFIG}" ]; then
    usage
fi

set -u

READS=$(echo "$@" | tr ' ' '\n')

SAMPLE="$(basename "$OUTPUTFOLDER")"

mkdir -p "$OUTPUTFOLDER"
Scont="$OUTPUTFOLDER/secontrol.txt"

echo 'Begin NewProject;' > "$Scont"
echo "File \"$OUTPUTFOLDER/${SAMPLE}.osprj\";" >> "$Scont"
echo "Options /Distributed=True;" >> "$Scont"
echo "End;" >> "$Scont"
echo "" >> "$Scont"

echo "Begin MapFusionReads /Namespace=NgsLib;" >> "$Scont"
echo "Files" >> "$Scont"
echo "\"$READS\";" >> "$Scont"


echo >> "$Scont"
cat "$CONFIG" | sed "s@OUTPUTFOLDER@$OUTPUTFOLDER@g" | sed "s@THREADS@$THREADS@g" >> "$Scont"
echo "Project $SAMPLE;" >> "$Scont"
echo "OutputFolder \"$OUTPUTFOLDER\";" >> "$Scont"
echo "End;" >> "$Scont"
echo "" >> "$Scont"

echo "Begin SaveProject;" >> "$Scont"
echo "Project $SAMPLE;" >> "$Scont"

echo "File \"$OUTPUTFOLDER/${SAMPLE}.osprj\";" >> "$Scont"
echo "End;" >> "$Scont"
echo "" >> "$Scont"
echo "Begin CloseProject;" >> "$Scont"
echo "Project $SAMPLE;" >> "$Scont"
echo "End;" >> "$Scont"


>&2 echo "project file completed:"
>&2 cat "$Scont"

cmd="mono /services/tools/oshell/10.0.1.29/oshell.exe --runscript FusionMapRef $OUTPUTFOLDER/secontrol.txt /scratch /bin/mono"
>&2 echo "starting oshell with command:"
>&2 echo "$cmd"
>&2 eval "$cmd"
