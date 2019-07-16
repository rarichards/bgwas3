import drmaa

def main():
    with drmaa.Session() as s:
        print('A DRMAA object was created')
        print('Supported contact strings: %s' % s.contact)
        print('Supported DRM systems: %s' % s.drmsInfo)
        print('Supported DRMAA implementations: %s' % s.drmaaImplementation)
        print('Exiting')

if __name__=='__main__':
    main()
