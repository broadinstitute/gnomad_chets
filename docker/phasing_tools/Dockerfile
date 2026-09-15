# Eagle, Beagle and SHAPEIT5 in one image, with GRCh38 genetic maps for all
# three, so the same region can be phased by each tool and scored the same way.
#
# One image rather than three keeps the comparison honest: same base OS, same
# bcftools, same maps derived from one source, so a difference between tools is
# a difference between tools.
#
# Build (no local Docker needed):
#   gcloud builds submit --project broad-mpg-gnomad \
#     --tag gcr.io/broad-mpg-gnomad/phasing_tools:5.1.1 docker/phasing_tools
FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl unzip \
        libgomp1 zlib1g \
    && rm -rf /var/lib/apt/lists/*

# SHAPEIT5 comes from bioconda: upstream's GitHub repository is disabled, so the
# release binaries are no longer downloadable from source. The conda package
# installs SHAPEIT5_phase_common, SHAPEIT5_phase_rare, SHAPEIT5_ligate and
# SHAPEIT5_switch. openjdk is for Beagle, which ships only as a jar.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        shapeit5=5.1.1 \
        bcftools=1.21 \
        htslib=1.21 \
        openjdk=21 \
    && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH

# Beagle 5.5. Pinned: the download is versioned in the filename and the
# unversioned name is not served, so this also documents exactly what ran.
RUN mkdir -p /opt/beagle \
    && curl -fsSL -o /opt/beagle/beagle.jar \
        https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar

# Eagle v2.4.1, matching the production local-ancestry pipeline. Its own
# genetic map ships inside the tarball.
RUN mkdir -p /opt/eagle \
    && curl -fsSL https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz \
        | tar -xz -C /opt/eagle --strip-components=1 \
    && chmod +x /opt/eagle/eagle

# Genetic maps. Beagle publishes a GRCh38 PLINK-format map; SHAPEIT5 wants a
# tab-separated "pos chr cM" file with a header. Deriving the SHAPEIT5 map from
# the Beagle one means both tools see the same recombination rates, so the
# comparison is not confounded by different maps.
#
# The chr_in_chrom_field variant matches GRCh38 contig naming (chr21, not 21),
# which is what the gnomAD VDS uses. Its files are named plink.chrchr21... --
# the doubled prefix is upstream's, not a typo here.
RUN mkdir -p /opt/maps/beagle /opt/maps/shapeit5 \
    && curl -fsSL -o /tmp/plink38.zip \
        https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip \
    && unzip -q /tmp/plink38.zip -d /tmp/plink38 \
    && for f in /tmp/plink38/chr_in_chrom_field/plink.chrchr*.GRCh38.map; do \
           c=$(basename "$f" | sed -E 's/^plink\.chr(chr[0-9XY]+)\.GRCh38\.map$/\1/'); \
           cp "$f" "/opt/maps/beagle/${c}.GRCh38.map"; \
           awk 'BEGIN{OFS="\t"; print "pos","chr","cM"} {print $4,$1,$3}' "$f" \
               > "/opt/maps/shapeit5/${c}.b38.gmap"; \
       done \
    && rm -rf /tmp/plink38.zip /tmp/plink38

ENV BEAGLE_JAR=/opt/beagle/beagle.jar \
    EAGLE_BIN=/opt/eagle/eagle \
    EAGLE_MAP=/opt/eagle/tables/genetic_map_hg38_withX.txt.gz \
    BEAGLE_MAP_DIR=/opt/maps/beagle \
    SHAPEIT_MAP_DIR=/opt/maps/shapeit5

# Fail the build rather than ship an image that is missing a tool. Each check
# names the tool it is on, so a failure identifies the culprit instead of
# collapsing to a bare exit code.
#
# Banners are captured to files rather than piped: the build shell runs with
# pipefail, and `tool | grep -q` kills the tool with SIGPIPE as soon as grep
# matches, which fails the build on a successful check. Eagle also rejects
# both --version and a bare invocation, so it is verified by its banner rather
# than by an exit code.
RUN set -eux \
    && SHAPEIT5_phase_common --help > /dev/null && echo "OK shapeit5 phase_common" \
    && SHAPEIT5_phase_rare --help > /dev/null && echo "OK shapeit5 phase_rare" \
    && SHAPEIT5_switch --help > /dev/null && echo "OK shapeit5 switch" \
    && SHAPEIT5_ligate --help > /dev/null && echo "OK shapeit5 ligate" \
    && bcftools --version > /tmp/v_bcftools.txt \
    && grep -q "^bcftools" /tmp/v_bcftools.txt && echo "OK bcftools" \
    && (java -jar "$BEAGLE_JAR" > /tmp/v_beagle.txt 2>&1 || true) \
    && grep -q "version 5.5" /tmp/v_beagle.txt && echo "OK beagle 5.5" \
    && ldd "$EAGLE_BIN" \
    && ("$EAGLE_BIN" > /tmp/v_eagle.txt 2>&1 || true) \
    && grep -q "Eagle v2.4.1" /tmp/v_eagle.txt && echo "OK eagle v2.4.1" \
    && test -s /opt/maps/beagle/chr21.GRCh38.map \
    && test -s /opt/maps/shapeit5/chr21.b38.gmap \
    && test -s "$EAGLE_MAP" \
    && head -1 /opt/maps/shapeit5/chr21.b38.gmap \
    && echo "OK maps" \
    && rm -f /tmp/v_*.txt

WORKDIR /work
