import logging
import xml.etree.ElementTree as ET
from datetime import datetime

class RunInfo:
    """Illumina sequencing run information.

    Attributes
    ----------
    run_id : str
        Illumina sequencing run ID.
    run_number : int
        Sequencing run number.
    flowcell_id : str
        Flowcell ID.
    instrument : str
        Illumina sequencing instrument ID.
    sequencing_date : datetime.date
        Date of the sequencing run.
    """
    def __init__(self, fname):
        self._parse(fname)
        self.logger = logging.getLogger(__name__)

    def _parse(self, fname):
        """Parse Illumina `RunInfo.xml`, set attributes of run information.
        """
        tree = ET.parse(fname)
        root = tree.getroot()

        run = root.find('Run')
        self.run_id = run.attrib['Id']
        self.run_number = int(run.attrib['Number'])
        self.flowcell_id = run.find('Flowcell').text
        self.instrument = run.find('Instrument').text
        self.sequencing_date = run.find('Date').text
        if len(self.sequencing_date) == 6:
            self.sequencing_date = datetime.strptime(self.sequencing_date, '%y%m%d').date()
        elif len(self.sequencing_date) == 8:
            self.sequencing_date = datetime.strptime(self.sequencing_date, '%Y%m%d').date()
        else:
            self.logger.warning(
                'Unrecognized sequencing date format: {}. Record raw string instead.'.format(self.sequencing_date)
            )
