from one.api import ONE
ONE.setup(base_url='https://openalyx.internationalbrainlab.org', silent=True)
one = ONE(password='international')
one = ONE()
eid = '4ecb5d24-f5cc-402c-be28-9d0f7cb14b3a'

passive_times = one.load_dataset(eid, '*passivePeriods*', collection='alf')
SP_times = passive_times['spontaneousActivity']