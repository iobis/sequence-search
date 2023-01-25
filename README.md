# sequence-search

## How to
### Generate bowtie2 database

First download all sequences and associated occurrences from the OBIS database to `data/sequences.csv` and `data/occurrences.csv`:

```sql
create temp table sequence_hashes as
select
	dna.id,
	dna.occurrence_id,
	MD5(dna.flat->>'DNA_sequence') as hash,
	dna.flat->>'DNA_sequence' as sequence
from dna;

create temp table sequence_occurrences as
select
	sh.hash,
	occurrence.decimallongitude as decimalLongitude,
	occurrence.decimallatitude as decimalLatitude,
	occurrence.dataset_id,
	aphia.classification->>'phylum' as phylum,
	aphia.classification->>'class' as class,
	aphia.classification->>'order' as order,
	aphia.classification->>'family' as family,
	aphia.classification->>'genus' as genus,
	aphia.record->>'scientificName' as scientificName,
	count(*)
from sequence_hashes sh
left join occurrence on sh.occurrence_id = occurrence.id
left join aphia on occurrence.aphia = aphia.id 
group by
	sh.hash,
	occurrence.decimallongitude,
	occurrence.decimallatitude,
	occurrence.dataset_id,
	aphia.classification->>'phylum',
	aphia.classification->>'class',
	aphia.classification->>'order',
	aphia.classification->>'family',
	aphia.classification->>'genus',
	aphia.record->>'scientificName';


select * from sequence_occurrences; -- occurrences.csv

select distinct on (hash) hash, sequence -- sequences.csv
from sequence_hashes
order by hash;
```

Build the fasta file, occurrence sqlite database, and bowtie2 database:

```sh
./build_db.sh
rsync -r data ***@***:/data/sequence-search/data
```

### Run the API

```sh
cd api
uvicorn main:app --reload --port 8086
# or
gunicorn --worker-class uvicorn.workers.UvicornWorker --bind '127.0.0.1:8086' --daemon main:app
```
