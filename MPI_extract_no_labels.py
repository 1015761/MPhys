import sys
import feature_extract as main
from utils import *
from os.path import exists
import pandas as pd
import numpy as np

try:
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	mpi_rank = comm.Get_rank()
	mpi_size = comm.Get_size()
	with_mpi = True
except ImportError:
	mpi_rank = 0
	mpi_size = 1
	with_mpi = False

mpi_root = 0

failed = []
times, star_radii, tic_list, labels, n_trans = get_tics(with_labels = False)
n_no_free = 0


if __name__ == '__main__':

	params = {}

	if mpi_rank == mpi_root:
		#Master node work

		nlc = 0
		for tic in times:
			for sec in times[tic]:
				nlc+=1

		if exists('/mnt/zfsusers/blakeland/Documents/MPhys/unseen.csv'): #Check this is right
			print("Existing manifest file found, will skip previously processed LCs and append to end of manifest file")
			sys.stdout.flush()
		else:
			print("Creating new manifest file")
			sys.stdout.flush()
			metadata_header = ['ID', 'TIC_ID','Sector', 'Times', 'Params', 'Uncs']
	
			with open('/mnt/zfsusers/blakeland/Documents/MPhys/unseen.csv', 'w') as f: # save in the photometry folder
				writer = csv.writer(f, delimiter=',')
				writer.writerow(metadata_header)


		if (not with_mpi) or (nlc == 1) or (mpi_size ==1):
			print("Not with MPI")
			print("Processing {} LCs".format(nlc))
			sys.stdout.flush()
			manifest_table = pd.read_csv('/mnt/zfsusers/blakeland/Documents/MPhys/unseen.csv')
			ID_exist = manifest_table['ID']
			for tic in times:
				for sec in times[tic]:
					#ID_sec = "s{}".format(10000 + int(sec))[1:]
					#ID_tic = str(1e16 + int(f))[1:]
					ID = '{}_{}'.format(tic, sec)


					if not np.isin(ID, ID_exist):
						print('Processing TIC {}'.format(tic))
						sys.stdout.flush()
						res = main.get_params(f, sec, 0)

						#if res['TIC'] < 0:
						#	print('Skipped TIC {}'.format(f))
						writer = csv.writer(file, delimiter=',')

						#ID_sec = "s{}".format(10000 + int(res['Sector']))[1:]
						#ID_tic = str(1e16 + int(res['TIC']))[1:]
						ID = '{}_{}'.format(res['TIC'], res['Sector'])


						#metadata_data.append(res['chunk'])
						writer.writerow(metadata_data)
						sys.stdout.flush()
						#	continue
						with open('/mnt/zfsusers/blakeland/Documents/MPhys/unseen.csv', 'a') as file: # save in the photometry folder
							writer = csv.writer(file, delimiter=',')

							#ID_sec = "s{}".format(10000 + int(res['Sector']))[1:]
							#ID_tic = str(1e16 + int(res['TIC']))[1:]
							ID = '{}_{}'.format(tic, sec)

							metadata_data = [ID]
							metadata_data.append(res['TIC_ID'])
							metadata_data.append(res['Sector'])
							metadata_data.append(res['Times'])
							metadata_data.append(res['Params'])
							metadata_data.append(res['Uncs'])


							#metadata_data.append(res['chunk'])
							writer.writerow(metadata_data)

					else:
						print('TIC {} already done - skipping'.format(f))
						sys.stdout.flush()

		else:
			#Master node delegating work to worker nodes


			sys.stdout.flush()

			if mpi_rank == 0:
				free_workers = list(range(1, mpi_size))
				active_workers = []
				n_finished = 0
				manifest_table = pd.read_csv('/mnt/zfsusers/blakeland/Documents/MPhys/unseen.csv')
				ID_exist = manifest_table['ID']

				while tic_list or active_workers:
					while tic_list and free_workers:
					#Send tics to workers
						tic = tic_list.pop()
						sec_list = list(times[tic])
						while sec_list and free_workers:
							sec = sec_list.pop()
							#ID_sec = "s{}".format(str(10000 + int(sec))[1:])
							#ID_tic = str(1000000000000000 + int(tic))[1:]
							ID = '{}_{}'.format(tic, sec)
							if np.isin(ID, ID_exist):
								print("TIC {}, sector {} light curve already analysed".format(tic, sec))
							else:
							
								if free_workers:
									w = free_workers.pop()
									comm.send(tic, dest = w, tag = 0)
									comm.send(sec, dest = w, tag = 1)
									active_workers.append(w)
									print('Worker {} is processing TIC {} in sector {}'.format(w, tic, sec))
									sys.stdout.flush()
								else:
									sec_list.append(sec)



								

					for w in active_workers:
						if comm.Iprobe(w, 2):
							res = comm.recv(source = w, tag = 2)
							print("Worker {} finished working".format(w))
							sys.stdout.flush()

							with open('/mnt/zfsusers/blakeland/Documents/MPhys/unseen.csv', 'a') as f:
								writer = csv.writer(f, delimiter=',')
								#if res['TIC'] != -99:
								#ID_sec = "s{}".format(10000 + int(res['Sector']))[1:]
								#ID_tic = str(1e16 + int(res['TIC']))[1:]
								ID = '{}_{}'.format(res['TIC'], res['Sector'])

								
								metadata_data = [ID]
								metadata_data.append(res['TIC_ID'])
								metadata_data.append(res['Sector'])
								metadata_data.append(res['Times'])
								metadata_data.append(res['Params'])

								try:
									metadata_data.append(res['Uncs'])
								except KeyError:
									print(res)

								#metadata_data.append(res['chunk'])
								writer.writerow(metadata_data)

							free_workers.append(w)
							active_workers.remove(w)
							n_finished += 1


				print('Completed {} of {} LCs'.format(n_finished, nlc))
				print('Had no free workers {} times'.format(n_no_free))
				sys.stdout.flush()


				for w in free_workers:
					comm.send(-1, dest = w, tag = 0)
				sys.stdout.flush()


	else:
		#Worker nodes
		while True:
			tic_worker = comm.recv(source = mpi_root, tag =0)
			if tic_worker == -1:
				break
			sec_worker = comm.recv(source = mpi_root, tag = 1)
			res = main.get_params(tic_worker, sec_worker, 0)
			comm.send(res, dest = mpi_root, tag = 2)



