# Troubleshooting commands
Here are several examples of code that does not work. Your job is to figure out why it doesn't work and correct the issue

```bash
ls /home/ubuntu/bam | grep ngm.rg.clean.bam > /home/ubuntu/tmp.bamfiles.list
java -jar $gatk \
   -T RealignTargetCreator \
   -R $ref \
   -I /home/ubuntu/${project}.bamfiles.list \
   -nt 2 \
   -log $log/tmp.RealignTargetCreator.log \
   -o /home/ubuntu/tmp.realign.intervals
```

