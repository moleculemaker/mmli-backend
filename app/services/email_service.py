import os
import smtplib
import jinja2
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.header import Header
from email.utils import formataddr

from config import app_config

SMTP_TIMEOUT_SECONDS = 30


class MailServerUnavailable(Exception):
    """The send failed before any recipient was involved: connect, greeting, or EHLO.

    Kept apart from every later failure because no recipient address can cause it, so
    callers may treat it as the server being down for everyone. Every later step can
    depend on the job's address, which is whatever was submitted with it: RCPT carries it
    directly, and so, indirectly, does MAIL FROM. Although MAIL FROM names only the
    configured sender, smtplib appends SIZE=<message length> when the server advertises
    SIZE, the message includes the To: header, and Job.email has no length limit. So an
    oversized address can get MAIL FROM refused.

    Known limitation: a server that answers EHLO and then stops responding still costs
    SMTP_TIMEOUT_SECONDS per send, since nothing here marks it unavailable. A timeout or
    disconnect at MAIL FROM (not a refusal) would be safe to treat as an outage, but that
    needs sendmail() split into mail()/rcpt()/data().
    """

class SingleEmailHeader(object):
    def __init__(self, recipients, email_params, template=''):
        self.recipients = recipients
        self.server = app_config['email']['server']
        self.from_email = app_config['email']['fromEmail']
        self.from_name = app_config['email']['fromName']
        ## Construct email o
        self.msg = MIMEMultipart('alternative')
        self.msg['Subject'] = email_params['subject']
        self.msg['From'] = formataddr((str(Header(self.from_name, 'utf-8')), self.from_email))
        self.msg['To'] = self.recipients
        ## Render email HTML content
        # self.html = self.render(os.path.join(os.path.dirname(__file__), 'templates', template), email_params)
        # self.msg.attach(MIMEText(self.html, 'html'))
        message_content_text = email_params['message']
        message_content = MIMEText(message_content_text, "plain")
        self.msg.attach(message_content)

    def render(self, tpl_path, email_params):
        path, filename = os.path.split(tpl_path)
        return jinja2.Environment(loader=jinja2.FileSystemLoader(path or './')).get_template(filename).render(email_params)

    def sendmail(self):
        # Without a timeout an unresponsive mail server blocks the caller indefinitely, and
        # the caller is usually the Kubernetes watcher thread.
        smtp_server = None
        try:
            smtp_server = smtplib.SMTP(self.server, port=25, timeout=SMTP_TIMEOUT_SECONDS)
            smtp_server.ehlo_or_helo_if_needed()
        except OSError as e:  # every smtplib exception is an OSError
            if smtp_server is not None:
                smtp_server.close()
            raise MailServerUnavailable(f'{self.server}: {e}') from e
        try:
            # The TO and CC header fields are populated by the header construction, and any additional recipient addresses are effectively BCC
            smtp_server.sendmail(self.from_email, self.recipients, self.msg.as_string())
        finally:
            # Not `with`: SMTP.__exit__ raises if QUIT is answered with anything but 221,
            # and by then the message has been delivered. Reporting that as a failed send
            # leaves the job due, so the watcher would send it again on every pass.
            try:
                smtp_server.quit()
            except (smtplib.SMTPException, OSError):
                smtp_server.close()

class EmailService:
    def __init__(self) -> None:
        pass

    def send_email(self, recipients=None, email_subject='', message_body='', template_file=None):
        if isinstance(recipients, list):
            recipients = recipients[0]
        email_params = {
            "subject": email_subject,
            "message": message_body,
        }
        email = SingleEmailHeader(recipients, email_params, template=template_file)
        # Raises on failure. Every caller catches it; the watcher depends on seeing it, since
        # a send reported as successful stamps job.notified_at and is never retried.
        email.sendmail()