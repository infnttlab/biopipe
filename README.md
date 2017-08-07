Ottimizzazione di una pipeline per l'analisi biofisica su cluster a basso consumo energetico

Prerequisiti:
- utilizzo della shell bash.


Installazione:
- eseguire l'installer install.sh presente nella directory;
- riavviare la bash.

Lo script comporta: 
- l'installazione di Miniconda3;
- l'installazione di alcune dipendenze necessarie al funzionamento del programma e all'analisi dei risultati(matplotlib, pandas, PyYAML, psutil e ipython);
- l'aggiunta di bioconda come canale per conda;
- la creazione di una cartella "data_ref" nella home directory contenente:
                                                                         - l'human reference: human_g1k_v37.fasta.gz
                                                                         - i due file fastq relativi al paziente analizzato:
                                                                                    - SRR1611178_1.fastq
                                                                                    - SRR1611178_2.fastq
La durata del processo dipende dal tipo di macchina, ma è mediamente compresa tra i 10 e i 15 minuti. I file contenuti nella cartella "data_ref" sono copiati dalla cartella "/mnt/avoton/biofisici/data_snakemake" .


Esecuzione del codice:
- eseguire script.sh presente nella directory.

Lo script è strutturato in tre fasi: 
- creazione dei subset da analizzare, tramite lo script split.py;
- analisi dei subset e creazione dei benchmarks per ogni regola eseguita, tramite l'esecuzione dello Snakefile;
- creazione di una tabella in formato csv che organizza i dati presenti nei benchmarks, tramite script_benchmark.py.

Come cambiare subsets:
- il programma split.py riceve due argomenti da linea di comando: il numero di read iniziale e il numero di read finale. 
Esempio) Un semplice subset potrebbe essere relativo alle prime 5000 reads del paziente e quindi generabile da: 
$  python ~/biopipe/split.py 0 5000 .
- La prima operazione di script.sh è la creazione dei subsets ed è gestita di default con un ciclo atto alla modifica dei parametri $start e $finish, che contrassegnano la prima read e l'ultima. Questa fase è soggetta a continue modifiche da parte dell'utente per poter analizzare subsets differenti.
- I subsets creati sono destinati ad una cartella nella directory biopipe denominata 'data/' che è rimossa al completamento dello script.

Analisi dei subsets:
- Per poter agevolare lo studio statistico dei risultati, lo Snakefile è eseguito di default 10 volte su ogni gruppo di subset, ovvero su coloro all'interno di 'data/'.

Risultato:
- l'oggetto di studio finale è la tabella ottenuta al termine dello script;
- La tabella è costituita dalle seguenti colonne:
[#cpu', 'Reference subject', 'cpu_type', 'h:m:s', 'io_in', 'io_out','max_pss', 'max_rss', 'max_uss', 'max_vms', 'mean_load', 'n_sim','rule', 's', 'threads']

Dove:
#cpu = numero di cpu utilizzate
Reference subject = specifica il subset
cpu_type = tipo di macchina utilizzata
n_sim = numero della simulazione
rule = specifica la regola
threads = numero di threads utilizzati
Le restanti colonne sono generate da snakemake nella creazione dei benchmarks e si affida alla documentazione ufficiale del programma la corretta descrizione.

Obiettivo:
- Ottenere conferma, dallo studio statistico dei dati organizzati nella tabella, che l'utilizzo di un cluster a basso consumo energico per l'analisi computazionale biofisica è potenzialmente più efficace e conveniente che l'uso dei moderni sistemi adoperati attualmente. 
