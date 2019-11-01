
# Practical differences Abel - Saga

- Combining storage: NIRD and computing power: SAGA
  - SAGA will be located in Trondheim. NIRD storage both in Tronheim and Tromsø (right?)
  - SAGA will be used for computing power only: at first there will be **NO BACKUP** from data hosted on Saga. It is planned to implement "slow file system where we could have some files at a later stage"

- **Disk space on Saga is much more limited on Saga than on Abel.** This means that we cannot have all our data permanently stored on Saga - only data required for ongoing analyses over some restricted period of time will have to be on Saga. 
> _in practice_: if you are finished with a certain type of data analysis (ie. your assembly), the the raw reads that you COPIED from NIRD to Saga to perform you analysis will have to be deleted/moved, and all your finished analyses will have to be moved to NIRD. If you do not do that, we might run out of space and nobody will be able to work

- **NIRD** data backup system is different than what was on Abel. On Abel, data were backed-up on tape. On NIRD this will be screenshot. 
- [ ] jeevan?
> if I understand well -> version control - Question: Does this have practical implications for users (said it will be more easy to messup)

- you will have the possibility to see your data directly from NIRD and Saga. (But this is not implemented yet). And you will **never be able to use data located on NIRD as input files**. 
  > This implies that you will have to transfer data back and forth from **NIRD** storage area to Saga to be able to do your analyses.
  
- **Disk areas on SAGA:**
  - `$HOME`
  - `$USERWORK`
  
  
- **Disk areas on NIRD:**
  
- [ ] jeevan, did I get this right? 
- There is a part on Saga that can be reserved for bioinformatic, part reserved for machine learning on SAGA (ie. if you need better floating point for your analyses)

# Administrators need to fix:
- find a way to have a common Conda environment 


# Where to find more detailed information

[SAGA documentation](https://documentation.sigma2.no/quick/saga.html)

[NIRD documentation](https://documentation.sigma2.no/storage/nird.html)

[Password changing](https://www.metacenter.no/user/login/)

[Abel to Saga course link](https://www.uio.no/english/services/it/research/events/Migrate-to-HPC.html) - includes HPC course.

Support: `support@metacenter.no` indicate which resource. See informations in documentation

# Summary: how to: main commands

## Login 

To saga 
```bash
ssh YOUR_USERNAME@saga.sigma2.no
```

- ?jeevan 
To NIRD - from within ABEL ? SAGA? 
```bash
login.nird.sigma2.no
```

## 
