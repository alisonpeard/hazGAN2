"""
SLURM fixer:

~/micromamba/envs/snakemake/

snakemake-executor-plugin-cluster-generic 1.0.9
snakemake-executor-plugin-slurm           1.1.0
snakemake-executor-plugin-slurm-jobstep   0.3.0
"""

# __init__.py/check_active_jobs
# sacct command
sacct_command = f"""sacct -X --parsable2 \
                --clusters all \
                --noheader --format=JobIdRaw,State \
                --starttime {sacct_starttime} \
                --endtime now --name {self.run_uuid}"""

# squeue alternative (1)
sacct_command = f"""squeue --noheader --all \
                --format '%i|%T' \
                --start {sacct_starttime} \
                --name {self.run_uuid}"""



# __init__.py/check_active_jobs
# -X: only show main job, no substeps
sacct_command = f"""sacct -X --parsable2 \
                --clusters all \
                --noheader --format=JobIdRaw,State \
                --starttime {sacct_starttime} \
                --endtime now --name {self.run_uuid}"""

# for better redability in verbose output
sacct_command = " ".join(shlex.split(sacct_command))

# this code is inspired by the snakemake profile:
# https://github.com/Snakemake-Profiles/slurm
for i in range(status_attempts):
    async with self.status_rate_limiter:
        (status_of_jobs, sacct_query_duration) = await self.job_stati(
            sacct_command
        )


# async def job_stati(self, command):
async def job_stati(self, command):
        """Obtain SLURM job status of all submitted jobs with sacct

        Keyword arguments:
        command -- a slurm command that returns one line for each job with:
                   "<raw/main_job_id>|<long_status_string>"
        """
        res = query_duration = None
        try:
            time_before_query = time.time()
            command_res = subprocess.check_output(
                command, text=True, shell=True, stderr=subprocess.PIPE
            )
            query_duration = time.time() - time_before_query
            self.logger.debug(
                f"The job status was queried with command: {command}\n"
                f"It took: {query_duration} seconds\n"
                f"The output is:\n'{command_res}'\n"
            )
            res = {
                # We split the second field in the output, as the State field
                # could contain info beyond the JOB STATE CODE according to:
                # https://slurm.schedmd.com/sacct.html#OPT_State
                entry[0]: entry[1].split(sep=None, maxsplit=1)[0]
                for entry in csv.reader(StringIO(command_res), delimiter="|")
            }
        except subprocess.CalledProcessError as e:
            error_message = e.stderr.strip()
            if "slurm_persist_conn_open_without_init" in error_message:
                self.logger.warning(
                    "The SLURM database might not be available ... "
                    f"Error message: '{error_message}'"
                    "This error message indicates that the SLURM database is currently "
                    "not available. This is not an error of the Snakemake plugin, "
                    "but some kind of server issue. "
                    "Please consult with your HPC provider."
                )
            else:
                self.logger.error(
                    f"The job status query failed with command '{command}'"
                    f"Error message: '{error_message}'"
                    "This error message is not expected, please report it back to us."
                )
            pass

        return (res, query_duration)

# sacct command
sacct_command = f"""sacct -X --parsable2 \
                --clusters all \
                --noheader --format=JobIdRaw,State \
                --starttime {sacct_starttime} \
                --endtime now --name {self.run_uuid}"""

# squeue alternative (1)
cmd = f"squeue --noheader --all --format '%i|%T' --start {sacct_starttime} --name {self.run_uuid}" \
