# first line: 391
@memory.cache
def run_mendel(data_dir, keep_results=True, **kwargs):
    #
    bootstrap_pedigrees = kwargs.pop('bootstrap_pedigrees', False)
    template_args = TemplateArgs(**kwargs)

    if not os.path.isdir(data_dir):
        data_dir = os.path.join(data_root_dir, data_dir)

    if not os.path.isdir(data_dir):
        raise ValueError(f'Non existing data dir {data_dir}')

    ped_name = os.path.basename(data_dir)
    cmd = generate_source_code(template_args)

    tdir = f'{data_dir}{hash(template_args)}'
    if bootstrap_pedigrees and os.path.isdir(tdir):
        # another process has processed this particular argument
        raise RuntimeError('Ignore bootstrap configuration.')

    print(f'Processing {ped_name} at {tdir}')
    if not os.path.isdir(tdir):
        os.makedirs(tdir)
        os.symlink(os.path.join(data_dir, 'locusbr.dat'), os.path.join(tdir, 'locusbr.dat'))
        # get number of pedigrees
        if bootstrap_pedigrees:
            with open(os.path.join(data_dir, 'pedigree.txt')) as ifile, \
                open(os.path.join(tdir, 'pedigree.txt'), 'w') as ofile:
                # first two lines
                ofile.write(ifile.readline())
                ofile.write(ifile.readline())
                pedigrees = []
                for line in ifile:
                    if len(line.split()) == 2:
                        # try to avoid duplicated name
                        nfam = line.split()[0]
                        pedigrees.append([f'{int(nfam): <5}PED{len(pedigrees)+ 1}\n'])
                        continue
                    pedigrees[-1].append(line)
                # write pedigrees.
                choices = np.random.choice(len(pedigrees), size=len(pedigrees), replace=True)
                choices.sort()
                for choice in choices:
                    ofile.write(''.join(pedigrees[choice]))
        else:
            os.symlink(os.path.join(data_dir, 'pedigree.txt'), os.path.join(tdir, 'pedigree.txt'))

    result_file = f'{tdir}/single_out.dat'
    proc = pexpect.spawn([cmd], cwd=tdir, timeout=500)

    # DO YOU WISH TO CHANGE TO BATCH MODE? [YES/NO]
    proc.expect(':')
    proc.sendline('no')
    # CHOOSE AN ITEM [1,...,21]:
    proc.expect(':')
    proc.sendline('21')
    # ANOTHER PROBLEM [YES/NO]:
    proc.expect(':')
    proc.sendline('no')
    # proc.interact()
    proc.expect(pexpect.EOF)

    # os.remove(cmd[:-2])
    # os.remove(cmd)
    res = extract_results(result_file) | asdict(template_args)
    # print(f'random samples kept in {tdir}')
    if not keep_results:
        shutil.rmtree(tdir)
    return res
