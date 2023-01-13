# sequence-search

## How to
### Setup

First download all sequences from the OBIS database to `data/sequences.csv`:

```sql
create temp table temp_sequences as
select
	'seq' || row_number() over () as id,
	dna.flat->>'target_gene' as target_gene,
	dna.flat->>'DNA_sequence' as sequence,
	string_agg(occurrence.id::text, '|') as occurrence_ids,
	string_agg(occurrence.aphia::text, '|') as aphia_ids
from dna
left join occurrence on dna.occurrence_id = occurrence.id
group by dna.flat->>'target_gene', dna.flat->>'DNA_sequence';
```

Build the fasta file, occurrence sqlite database, and bowtie2 database:

```sh
python generate_files.py
./build_db.sh
```

### Run the API

```sh
cd api
uvicorn main:app --reload --port 8086
# or
gunicorn --worker-class uvicorn.workers.UvicornWorker --bind '127.0.0.1:8086' --daemon main:app
```
