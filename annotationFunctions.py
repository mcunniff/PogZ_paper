
def nonOptoZoneAnnotations(filename, sheet, zoneNames, zoneTags, anymazestart = 0, output = 'zones.txt',lengthFilter = True):
    """Annotations based on position of animal in number of zones, no light information"""
    
    df = pd.read_excel(filename,sheet)
    records = []
    allTags = []
    for i in range(0,len(zoneNames)):
        startcol = zoneNames[i] + 'Start'
        stopcol = zoneNames[i] + 'Stop'

        # finds all stop & start times from anymaze position file
        starts = df.iloc[df[startcol].nonzero()][startcol]
        starts = starts.dropna().as_matrix()
        stops = df.iloc[df[stopcol].nonzero()][stopcol]
        stops = stops.dropna().as_matrix()

        # checks for trailing start event, eliminates if present 
        if starts[-1] > stops[-1]:
            starts = starts[0:(len(starts)-1)]

        # checks if both start and stop are in a given epoch, assigns to that epoch if so
        # aside from trailing start, doesn't check for proper matching at this time
        events = pd.DataFrame({'Start':starts, 'Stop':stops})
        events['Zone'] = zoneTags[i]

        records.append(events)
        
    cols = ['Start','Stop','Zone']
    #combine all individual zones into one master record
    allRec = pd.concat(records, ignore_index = True)

	# # eliminates events that are <2.5 seconds because of FFT window
    if lengthFilter == True:
        allRec['Length'] = allRec['Stop'] - allRec['Start']
        allRec = allRec[(allRec['Length'] > 2.5)]

    # allRec = allRec.sort_values(['Start'])[cols].reset_index(drop = True)
    allRec['Start'] = allRec['Start'] + anymazestart
    allRec['Stop'] = allRec['Stop'] + anymazestart
    allRec.to_csv(output, index = False, header = False, columns = ['Start','Stop','Zone'],sep=' ')



def nonOptoZoneFileList(positionFile, positionSheets, outputs, data, finalOutput, zoneNames, zoneTags, anymaze_starts,lengthFilter = True):
    """generates file list for matlab input using classification by single zone (in or out)"""
    try:
        len(anymaze_starts)
    except NameError:
        anymaze_starts = zeros(len(positionSheets))
    for k in range(0,len(positionSheets)):
        print positionSheets[k]
        nonOptoZoneAnnotations(positionFile,positionSheets[k], zoneNames, zoneTags, anymaze_starts[k], outputs[k],lengthFilter)
    cols = ['Files','Annotations']
    records = pd.DataFrame({'Files':data,'Annotations':outputs})
    records = records.sort_values(['Files'])[cols].reset_index(drop = True)
    records.to_csv(finalOutput, sep='\t', index=False, header=False)


def annotateRunningTimes(filename, sheet, colName, anymazestart = 0, output = 'zones.txt'):
    """Selects times around transitions as proxy for time animal is moving"""
    import pandas as pd
    df = pd.read_excel(filename,sheet)
    records = []
    allTags = []
    startcol = colName + 'Start'
    stopcol = colName + 'Stop'

    # finds all stop & start times from anymaze position file, gets indicies to compare 
    starts = df.iloc[df[startcol].nonzero()][startcol]
    starts = starts.dropna()
    stops = df.iloc[df[stopcol].nonzero()][stopcol]
    stops = stops.dropna()

    starts = starts.values
    stops = stops.values

    if starts[-1] > stops[-1]:
        starts = starts[0:(len(starts)-1)]

    events = pd.DataFrame({'Start':starts,'Stop':stops})
    events.drop(events[(events['Start'] == events['Stop'])].index, inplace=True)
    events['Start'] = events['Start'] + anymazestart

    runStart = events['Start'] - 1.5;
    runStop = events['Start'] + 1.5;


    runEvents = pd.DataFrame({'Start': runStart, 'Stop': runStop})
    runEvents['Type'] = 0;

    runEvents.drop(runEvents[runEvents['Start'] < 0].index, inplace=True)

    allRec = runEvents[['Start','Stop','Type']]
    allRec.to_csv(output, index = False, header = False, columns = ['Start','Stop','Type'],sep=' ')

def generateRunningTimeFileList(positionFile, positionSheets, outputs, data, finalOutput, colName, anymaze_starts):
    import pandas as pd
    """takes list of input files, names for annotation files, generates annotations, & outputs filelist"""
    for k in range(0,len(data)):
        print(positionSheets[k])
        annotateRunningTimes(positionFile, positionSheets[k], colName, anymaze_starts[k], outputs[k])    
    cols = ['Files','Annotations']
    records = pd.DataFrame({'Files':data,'Annotations':outputs})
    records = records.sort_values(['Files'])[cols].reset_index(drop = True)
    records.to_csv(finalOutput, sep=' ', index=False, header=False)

def pogzRunTypes(filename, sheet, colName, compareCol, anymazestart = 0, output = 'zones.txt'):
    """Classifies entries to zone based on whether additional zone is entered (eg stim zone & open arms)"""
    # 0: not in zone, light off
    # 1: not in zone, light on
    # 2: in zone, light off
    # 3: in zone, light on
        # finds all stop & start times from anymaze position file, gets indicies to compare 
    starts = df.iloc[df[startcol].nonzero()][startcol]
    starts = starts.dropna()
    stops = df.iloc[df[stopcol].nonzero()][stopcol]
    stops = stops.dropna()
    
    startIndicies = starts.index.values 
    stopIndicies = stops.index.values
    starts = starts.values
    stops = stops.values

    if starts[-1] > stops[-1]:
        starts = starts[0:(len(starts)-1)]
        startIndicies = startIndicies[0:(len(starts)-1)]

    # for each entry, sums up all events in compare column (should be "In Open Arm" or similar one/zero array)
    eventSums = []
    for n in range(0,len(startIndicies)):
        openSum = int(df[compareCol].loc[startIndicies[n]:stopIndicies[n]].sum())
        eventSums.append(openSum)

    # aligns starts & stops and classifies based on eventSum value
    events = pd.DataFrame({'Start':starts,'Stop':stops,'Compare': eventSums})
    events['Run Type'] = np.where(events['Compare'] > 2,1,0)

    events['Type'] = events['Run Type']
    events['Start'] = events['Start'] + anymazestart
    events['Stop'] = events['Stop'] + anymazestart

    allRec = events[['Start','Stop','Type']]
    allRec.to_csv(output, index = False, header = False, columns = ['Start','Stop','Type'],sep=' ')

def runTypeFileList(positionFile, positionSheets, outputs, data, finalOutput, colName, compareCol, anymaze_starts):
    """takes list of input files, names for annotation files, generates annotations, & outputs filelist"""
    for k in range(0,len(data)):
        print(positionSheets[k])
        pogzRunTypes(positionFile, positionSheets[k], colName, compareCol, anymaze_starts[k], outputs[k])    
    cols = ['Files','Annotations']
    records = pd.DataFrame({'Files':data,'Annotations':outputs})
    records = records.sort_values(['Files'])[cols].reset_index(drop = True)
    records.to_csv(finalOutput, sep=' ', index=False, header=False)