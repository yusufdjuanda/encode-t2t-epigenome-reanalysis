# Resources

Keep small manifests and checksums here. Do not commit reference FASTA files,
Bowtie2 indexes, or large annotation tracks.

For every external resource, record:

- source URL;
- accession or release identifier;
- download date;
- SHA-256 checksum;
- coordinate system and reference assembly;
- license or usage restrictions.

The paper-exact analysis must freeze the same T2T-CHM13 release and associated
annotations used by the publication. A modern analysis should use a separate
resource directory and must never silently replace the paper-exact reference.
