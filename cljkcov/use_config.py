import yaml


def load_yaml_config(ymlfile, print_sim_info=True) :
    ''' Loads the yml file for simulation and jackknife configuration'''

    with open(ymlfile,'r') as f:
        data = yaml.full_load(f)

    print('Using config for: ', data.get('simulation_name'))
    print('Info: ', data.get('simulation_description'))


    assert "inputfile" in data, "No inputfile defined in yml."


    return data
