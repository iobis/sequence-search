# sequence-search

## How to

First download all sequences from the OBIS database:

```sql
create temp table temp_sequences as
select
	'seq' || row_number() over () as id,
	dna.flat->>'target_gene' as target_gene,
	dna.flat->>'DNA_sequence' as sequence,
	string_agg(id::text, '|') as ids
from dna
group by dna.flat->>'target_gene', dna.flat->>'DNA_sequence';
```

Build the fasta file:

```sh
python generate_fasta.py
```

Build the bowtie2 database:

```sh
./build_db.sh
```

Generate an input fasta file:

```sh
python generate_input.py
```

Test search:

```sh
./search.sh
```
